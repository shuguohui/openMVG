// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/tracks/tracks.hpp"
#include "openMVG/types.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

#include <ceres/types.h>

#include <iostream>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG{
namespace sfm{

using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;


// Sorts the grid cell elements by the track statistics, which will sort first
// by the (truncated) track length, then by the mean reprojection error.
bool CompareGridCellElements(
	const std::pair<IndexT, GlobalSfMReconstructionEngine_RelativeMotions::TrackStatistics>& element1,
	const std::pair<IndexT, GlobalSfMReconstructionEngine_RelativeMotions::TrackStatistics>& element2) {
	return element1.second < element2.second;
}

GlobalSfMReconstructionEngine_RelativeMotions::GlobalSfMReconstructionEngine_RelativeMotions(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory), sLogging_file_(sloggingFile)
{

  if (!sLogging_file_.empty())
  {
    // setup HTML logger
    html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("GlobalReconstructionEngine SFM report.");
    html_doc_stream_->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("GlobalSfMReconstructionEngine_RelativeMotions")));
    html_doc_stream_->pushInfo("<hr>");

    html_doc_stream_->pushInfo( "Dataset info:");
    html_doc_stream_->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }

  // Set default motion Averaging methods
  eRotation_averaging_method_ = ROTATION_AVERAGING_L2;
  eTranslation_averaging_method_ = TRANSLATION_AVERAGING_L1;
}

GlobalSfMReconstructionEngine_RelativeMotions::~GlobalSfMReconstructionEngine_RelativeMotions()
{
  if (!sLogging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(sLogging_file_.c_str());
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetFeaturesProvider(Features_Provider * provider)
{
  features_provider_ = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetMatchesProvider(Matches_Provider * provider)
{
  matches_provider_ = provider;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetRotationAveragingMethod
(
  ERotationAveragingMethod eRotationAveragingMethod
)
{
  eRotation_averaging_method_ = eRotationAveragingMethod;
}

void GlobalSfMReconstructionEngine_RelativeMotions::SetTranslationAveragingMethod
(
  ETranslationAveragingMethod eTranslationAveragingMethod
)
{
  eTranslation_averaging_method_ = eTranslationAveragingMethod;
}

bool GlobalSfMReconstructionEngine_RelativeMotions::Process() {


	std::string temp_sfm_data = stlplus::create_filespec(sOut_directory_, "temp_sfm_data.bin");
	if (b_save_intermediate_result_ && stlplus::file_exists(temp_sfm_data))
	{
		Load(sfm_data_, temp_sfm_data, ALL);
	}
	else
	{
		//-------------------
		// Keep only the largest biedge connected subgraph
		//-------------------
		{
			const Pair_Set pairs = matches_provider_->getPairs();
			const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
			if (set_remainingIds.empty())
			{
				std::cout << "Invalid input image graph for global SfM" << std::endl;
				return false;
			}
			KeepOnlyReferencedElement(set_remainingIds, matches_provider_->pairWise_matches_);
		}

		openMVG::rotation_averaging::RelativeRotations relatives_R;
		Compute_Relative_Rotations(relatives_R);

		Hash_Map<IndexT, Mat3> global_rotations;
		if (!Compute_Global_Rotations(relatives_R, global_rotations))
		{
			std::cerr << "GlobalSfM:: Rotation Averaging failure!" << std::endl;
			return false;
		}

		matching::PairWiseMatches  tripletWise_matches;
		if (!Compute_Global_Translations(global_rotations, tripletWise_matches))
		{
			std::cerr << "GlobalSfM:: Translation Averaging failure!" << std::endl;
			return false;
		}
		if (!Compute_Initial_Structure(tripletWise_matches))
		{
			std::cerr << "GlobalSfM:: Cannot initialize an initial structure!" << std::endl;
			return false;
		}
		if (!SelectGoodTracksForBundleAdjustment())
		{
			std::cerr << "GlobalSfM:: Select GoodTracks For BundleAdjustment failure!" << std::endl;
			return false;
		}
		if(b_save_intermediate_result_)
			Save(sfm_data_, temp_sfm_data, ALL);
	}



	if (!Adjust())
	{
		std::cerr << "GlobalSfM:: Non-linear adjustment failure!" << std::endl;
		return false;
	}

	//-- Export statistics about the SfM process
	if (!sLogging_file_.empty())
	{
		using namespace htmlDocument;
		std::ostringstream os;
		os << "Structure from Motion statistics.";
		html_doc_stream_->pushInfo("<hr>");
		html_doc_stream_->pushInfo(htmlMarkup("h1", os.str()));

		os.str("");
		os << "-------------------------------" << "<br>"
			<< "-- View count: " << sfm_data_.GetViews().size() << "<br>"
			<< "-- Intrinsic count: " << sfm_data_.GetIntrinsics().size() << "<br>"
			<< "-- Pose count: " << sfm_data_.GetPoses().size() << "<br>"
			<< "-- Track count: " << sfm_data_.GetLandmarks().size() << "<br>"
			<< "-------------------------------" << "<br>";
		html_doc_stream_->pushInfo(os.str());
	}

	return true;
}

/// Compute from relative rotations the global rotations of the camera poses
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Rotations
(
  const rotation_averaging::RelativeRotations & relatives_R,
  Hash_Map<IndexT, Mat3> & global_rotations
)
{
  if(relatives_R.empty())
    return false;
  // Log statistics about the relative rotation graph
  {
    std::set<IndexT> set_pose_ids;
    for (const auto & relative_R : relatives_R)
    {
      set_pose_ids.insert(relative_R.i);
      set_pose_ids.insert(relative_R.j);
    }

    std::cout << "\n-------------------------------" << "\n"
      << " Global rotations computation: " << "\n"
      << "  #relative rotations: " << relatives_R.size() << "\n"
      << "  #global rotations: " << set_pose_ids.size() << std::endl;
  }

  // Global Rotation solver:
  const ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod =
    TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR;
    //TRIPLET_ROTATION_INFERENCE_NONE;

  system::Timer t;
  GlobalSfM_Rotation_AveragingSolver rotation_averaging_solver;
  const bool b_rotation_averaging = rotation_averaging_solver.Run(
    eRotation_averaging_method_, eRelativeRotationInferenceMethod,
    relatives_R, global_rotations);

  std::cout
    << "Found #global_rotations: " << global_rotations.size() << "\n"
    << "Timing: " << t.elapsed() << " seconds" << std::endl;


  if (b_rotation_averaging)
  {
    // Compute & display rotation fitting residual errors
    std::vector<float> vec_rotation_fitting_error;
    vec_rotation_fitting_error.reserve(relatives_R.size());
    for (const auto & relative_R : relatives_R)
    {
      const Mat3 & Rij = relative_R.Rij;
      const IndexT i = relative_R.i;
      const IndexT j = relative_R.j;
      if (global_rotations.count(i)==0 || global_rotations.count(j)==0)
        continue;
      const Mat3 & Ri = global_rotations[i];
      const Mat3 & Rj = global_rotations[j];
      const Mat3 eRij(Rj.transpose()*Rij*Ri);
      const double angularErrorDegree = R2D(getRotationMagnitude(eRij));
      vec_rotation_fitting_error.push_back(angularErrorDegree);
    }

    if (!vec_rotation_fitting_error.empty())
    {
      const float error_max = *max_element(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
      Histogram<float> histo(0.0f,error_max, 20);
      histo.Add(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
      std::cout
        << "\nRelative/Global degree rotations residual errors {0," << error_max<< "}:"
        << histo.ToString() << std::endl;
      {
        Histogram<float> histo(0.0f, 5.0f, 20);
        histo.Add(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
        std::cout
          << "\nRelative/Global degree rotations residual errors {0,5}:"
          << histo.ToString() << std::endl;
      }
      std::cout << "\nStatistics about global rotation evaluation:" << std::endl;
      minMaxMeanMedian<float>(vec_rotation_fitting_error.begin(), vec_rotation_fitting_error.end());
    }

    // Log input graph to the HTML report
    if (!sLogging_file_.empty() && !sOut_directory_.empty())
    {
      // Log a relative pose graph
      {
        std::set<IndexT> set_pose_ids;
        Pair_Set relative_pose_pairs;
        for (const auto & view : sfm_data_.GetViews())
        {
          const IndexT pose_id = view.second->id_pose;
          set_pose_ids.insert(pose_id);
        }
        const std::string sGraph_name = "global_relative_rotation_pose_graph_final";
        graph::indexedGraph putativeGraph(set_pose_ids, rotation_averaging_solver.GetUsedPairs());
        graph::exportToGraphvizData(
          stlplus::create_filespec(sOut_directory_, sGraph_name),
          putativeGraph);

        using namespace htmlDocument;
        std::ostringstream os;

        os << "<br>" << sGraph_name << "<br>"
           << "<img src=\""
           << stlplus::create_filespec(sOut_directory_, sGraph_name, "svg")
           << "\" height=\"600\">\n";

        html_doc_stream_->pushInfo(os.str());
      }
    }
  }
  return b_rotation_averaging;
}

/// Compute/refine relative translations and compute global translations
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Global_Translations
(
  const Hash_Map<IndexT, Mat3> & global_rotations,
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Translation averaging (compute translations & update them to a global common coordinates system)
  GlobalSfM_Translation_AveragingSolver translation_averaging_solver;
  const bool bTranslationAveraging = translation_averaging_solver.Run(
    eTranslation_averaging_method_,
    sfm_data_,
    features_provider_,
    matches_provider_,
    global_rotations,
    tripletWise_matches);

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "cameraPath_translation_averaging", "ply"),
      ESfM_Data(EXTRINSICS));
  }

  return bTranslationAveraging;
}

/// Compute the initial structure of the scene
bool GlobalSfMReconstructionEngine_RelativeMotions::Compute_Initial_Structure
(
  matching::PairWiseMatches & tripletWise_matches
)
{
  // Build tracks from selected triplets (Union of all the validated triplet tracks (_tripletWise_matches))
  {
    using namespace openMVG::tracks;
    TracksBuilder tracksBuilder;
#if defined USE_ALL_VALID_MATCHES // not used by default
    matching::PairWiseMatches pose_supported_matches;
    for (const std::pair< Pair, IndMatches > & match_info :  matches_provider_->pairWise_matches_)
    {
      const View * vI = sfm_data_.GetViews().at(match_info.first.first).get();
      const View * vJ = sfm_data_.GetViews().at(match_info.first.second).get();
      if (sfm_data_.IsPoseAndIntrinsicDefined(vI) && sfm_data_.IsPoseAndIntrinsicDefined(vJ))
      {
        pose_supported_matches.insert(match_info);
      }
    }
    tracksBuilder.Build(pose_supported_matches);
#else
    // Use triplet validated matches
    tracksBuilder.Build(tripletWise_matches);
#endif
    tracksBuilder.Filter(3);
    STLMAPTracks map_selectedTracks; // reconstructed track (visibility per 3D point)
    tracksBuilder.ExportToSTL(map_selectedTracks);

    // Fill sfm_data with the computed tracks (no 3D yet)
    Landmarks & structure = sfm_data_.structure;
    IndexT idx(0);
    for (STLMAPTracks::const_iterator itTracks = map_selectedTracks.begin();
      itTracks != map_selectedTracks.end();
      ++itTracks, ++idx)
    {
      const submapTrack & track = itTracks->second;
      structure[idx] = Landmark();
      Observations & obs = structure.at(idx).obs;
      for (submapTrack::const_iterator it = track.begin(); it != track.end(); ++it)
      {
        const size_t imaIndex = it->first;
        const size_t featIndex = it->second;
        const PointFeature & pt = features_provider_->feats_per_view.at(imaIndex)[featIndex];
        obs[imaIndex] = Observation(pt.coords().cast<double>(), featIndex);
      }
    }

    std::cout << std::endl << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats:
      //    - number of images
      //    - number of tracks
      std::set<uint32_t> set_imagesId;
      TracksUtilsMap::ImageIdInTracks(map_selectedTracks, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<uint32_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<uint32_t, uint32_t> map_Occurence_TrackLength;
      TracksUtilsMap::TracksLength(map_selectedTracks, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (const auto & iter : map_Occurence_TrackLength)  {
        osTrack << "\t" << iter.first << "\t" << iter.second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }

  // Compute 3D position of the landmark of the structure by triangulation of the observations
  {
    openMVG::system::Timer timer;

    const IndexT trackCountBefore = sfm_data_.GetLandmarks().size();
    SfM_Data_Structure_Computation_Blind structure_estimator(true);
    structure_estimator.triangulate(sfm_data_);

    std::cout << "\n#removed tracks (invalid triangulation): " <<
      trackCountBefore - IndexT(sfm_data_.GetLandmarks().size()) << std::endl;
    std::cout << std::endl << "  Triangulation took (s): " << timer.elapsed() << std::endl;

    // Export initial structure
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "initial_structure", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }
  return !sfm_data_.structure.empty();
}

void GlobalSfMReconstructionEngine_RelativeMotions::ComputeTrackStatistics(const int long_track_length_threshold,
	std::unordered_map<IndexT, TrackStatistics>& track_statistics)
{
	const Landmarks& landmarks = sfm_data_.structure;
	track_statistics.reserve(landmarks.size());

	// Compute the track statistics for each track.
	for (Landmarks::const_iterator citer =  landmarks.begin();
		citer != landmarks.end();
		++citer)
	{
		const IndexT trackID = citer->first;
		const Landmark& landmark = citer->second;
		double sq_reprojection_error_sum = 0.0;
		int num_valid_reprojections = 0;
		// Compute the sq reprojection error for each view that observes the track
		// and it it to the accumulating sum.
		for (Observations::const_iterator citer_obs = landmark.obs.begin();
			citer_obs != landmark.obs.end();
			++citer_obs)
		{
			const IndexT viewID = citer_obs->first;
			View* pView = sfm_data_.views[viewID].get();
			cameras::IntrinsicBase* pIntrinsic = sfm_data_.intrinsics[pView->id_intrinsic].get();
			Vec2 p = pIntrinsic->project(sfm_data_.poses[viewID], landmark.X);
			sq_reprojection_error_sum += ( p - citer_obs->second.x).squaredNorm();
		}
		

		// Set the track statistics in the output map.
		const int truncated_track_length =
			std::min((int)landmark.obs.size(), long_track_length_threshold);
		const double mean_sq_reprojection_error =
			sq_reprojection_error_sum / static_cast<double>(landmark.obs.size());
		track_statistics.emplace(trackID,std::make_pair(truncated_track_length, mean_sq_reprojection_error));
	}
}

void GlobalSfMReconstructionEngine_RelativeMotions::SelectBestTracksFromEachImageGridCell(
	IndexT viewID,
	int image_grid_cell_size,
	const std::unordered_set<IndexT>& trackids,
	const std::unordered_map<IndexT, TrackStatistics>& track_statistics,
	std::unordered_set<IndexT>& tracks_to_optimize)
{
	const double inv_grid_cell_size = 1.0 / image_grid_cell_size;

	// Hash each feature into a grid cell.
	ImageGrid image_grid;
	
	for (const IndexT track_id : trackids)
	{
		const Landmark& landmark = sfm_data_.structure[track_id];
		const Vec2& pos = landmark.obs.at(viewID).x;
		const TrackStatistics& current_track_statistics = track_statistics.at(track_id);
		int64_t dx = pos[0] * inv_grid_cell_size;
		int64_t dy = pos[1] * inv_grid_cell_size;
		int64_t gridID = (dx << 32) | dy;

		image_grid[gridID].emplace_back(track_id, current_track_statistics);
	}

	// Select the best feature from each grid cell and add it to the tracks to
	// optimize.
	for (auto& grid_cell : image_grid) 
	{
		// Order the features in each cell by track length first, then mean
		// reprojection error.
		const GridCellElement& grid_cell_element =
			*std::min_element(grid_cell.second.begin(),
				grid_cell.second.end(),
				CompareGridCellElements);

		// Insert the track id in to the tracks to optimize.
		tracks_to_optimize.emplace(grid_cell_element.first);
	}
}
void GlobalSfMReconstructionEngine_RelativeMotions::SelectTopRankedTracksInView(
	const std::unordered_map<IndexT, TrackStatistics>& track_statistics,
	const std::unordered_set<IndexT>& trackids,
	IndexT viewID, int min_num_optimized_tracks_per_view,
	std::unordered_set<IndexT>& tracks_to_optimize)
{
	int num_optimized_tracks = 0;
	int num_estimated_tracks = 0;

	
	std::vector<GridCellElement> ranked_candidate_tracks;
	for (const IndexT track_id : trackids)
	{
		const Landmark& landmark = sfm_data_.structure[track_id];

		// We only reach this point if the track is estimated.
		++num_estimated_tracks;

		// If the track is already slated for optimization, increase the count of
		// optimized features.
		if (tracks_to_optimize.find(track_id) != tracks_to_optimize.end())
		{
			++num_optimized_tracks;
			// If the number of optimized_tracks is greater than the minimum then we
			// can return early since we know that no more features need to added for
			// this view.
			if (num_optimized_tracks >= min_num_optimized_tracks_per_view) 
			{
				return;
			}
		}
		else 
		{
			// If the track is not already set to be optimized then add it to the list
			// of candidate tracks.
			ranked_candidate_tracks.emplace_back(track_id, track_statistics.at(track_id));
		}
	}

	// We only reach this point if the number of optimized tracks is less than the
	// minimum. If that is the case then we add the top candidate features until
	// the minimum number of features observed is met.
	if (num_optimized_tracks != num_estimated_tracks) 
	{
		// Select how many tracks to add. If we need more tracks than are estimated
		// then we simply add all remaining features.
		const int num_optimized_tracks_needed =
			std::min(min_num_optimized_tracks_per_view - num_optimized_tracks,
				num_estimated_tracks - num_optimized_tracks);
		std::partial_sort(
			ranked_candidate_tracks.begin(),
			ranked_candidate_tracks.begin() + num_optimized_tracks_needed,
			ranked_candidate_tracks.end());
		// Add the candidate tracks to the list of tracks to be optimized.
		for (int i = 0; i < num_optimized_tracks_needed; i++) {
			tracks_to_optimize.emplace(ranked_candidate_tracks[i].first);
		}
	}
}

bool GlobalSfMReconstructionEngine_RelativeMotions::SelectGoodTracksForBundleAdjustment()
{
	const int long_track_length_threshold = 10;
	const int image_grid_cell_size = 100;
	const int min_num_optimized_tracks_per_view = 0;
	// Compute the track mean reprojection errors.
	std::unordered_map<IndexT, TrackStatistics> track_statistics;
	ComputeTrackStatistics(long_track_length_threshold,track_statistics);


	//collect each view's tracks
	std::unordered_map<IndexT, std::unordered_set<IndexT> > view_tracks_map;
	for (Landmarks::const_iterator citer = sfm_data_.structure.begin();
		citer != sfm_data_.structure.end();
		++citer)
	{
		const IndexT trackID = citer->first;
		const Landmark& landmark = citer->second;
		for (Observations::const_iterator citer_obs = landmark.obs.begin();
			citer_obs != landmark.obs.end();
			++citer_obs)
		{
			const IndexT viewID = citer_obs->first;
			view_tracks_map[viewID].insert(trackID);
		}
	}
	// For each image, divide the image into a grid and choose the highest quality
	// tracks from each grid cell. This encourages good spatial coverage of tracks
	// within each image.
	std::unordered_set<IndexT> tracks_to_optimize;
	for(Views::const_iterator citer = sfm_data_.views.begin();
		citer != sfm_data_.views.end();
		++citer)
	{
		const IndexT viewID = citer->first;
		SelectBestTracksFromEachImageGridCell(viewID, image_grid_cell_size,view_tracks_map[viewID], track_statistics, tracks_to_optimize);
	}


	// To this point, we have only added features that have as full spatial
	// coverage as possible within each image but we have not ensured that each
	// image is constrainted by at least K features. So, we cycle through all
	// views again and add the top M tracks that have not already been added.
	for (Views::const_iterator citer = sfm_data_.views.begin();
		citer != sfm_data_.views.end();
		++citer)
	{
		const IndexT viewID = citer->first;

		// If this view is not constrained by enough optimized tracks, add the top
		// ranked features until there are enough tracks constraining the view.
		SelectTopRankedTracksInView(track_statistics, view_tracks_map[viewID],viewID,min_num_optimized_tracks_per_view,tracks_to_optimize);
	}

	Landmarks newLandmarks;
	for (auto track_id :  tracks_to_optimize)
	{
		newLandmarks[track_id] = sfm_data_.structure[track_id];
	}
	sfm_data_.structure.swap(newLandmarks);

	return true;
}
// Adjust the scene (& remove outliers)
bool GlobalSfMReconstructionEngine_RelativeMotions::Adjust()
{
  // Refine sfm_scene (in a 3 iteration process (free the parameters regarding their incertainty order)):

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  // - refine only Structure and translations
  bool b_BA_Status = bundle_adjustment_obj.Adjust
    (
      sfm_data_,
      Optimize_Options(
        Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
        Extrinsic_Parameter_Type::ADJUST_TRANSLATION, // Rotations are held as constant
        Structure_Parameter_Type::ADJUST_ALL,
        Control_Point_Parameter(),
        this->b_use_motion_prior_)
    );
  if (b_BA_Status)
  {
    if (!sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_00_refine_T_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }

    // - refine only Structure and Rotations & translations
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          Intrinsic_Parameter_Type::NONE, // Intrinsics are held as constant
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_01_refine_RT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  if (b_BA_Status && ReconstructionEngine::intrinsic_refinement_options_ != Intrinsic_Parameter_Type::NONE) {
    // - refine all: Structure, motion:{rotations, translations} and optics:{intrinsics}
    b_BA_Status = bundle_adjustment_obj.Adjust
      (
        sfm_data_,
        Optimize_Options(
          ReconstructionEngine::intrinsic_refinement_options_,
          Extrinsic_Parameter_Type::ADJUST_ALL,
          Structure_Parameter_Type::ADJUST_ALL,
          Control_Point_Parameter(),
          this->b_use_motion_prior_)
      );
    if (b_BA_Status && !sLogging_file_.empty())
    {
      Save(sfm_data_,
        stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_02_refine_KRT_Xi", "ply"),
        ESfM_Data(EXTRINSICS | STRUCTURE));
    }
  }

  // Remove outliers (max_angle, residual error)
  const size_t pointcount_initial = sfm_data_.structure.size();
  RemoveOutliers_PixelResidualError(sfm_data_, 4.0);
  const size_t pointcount_pixelresidual_filter = sfm_data_.structure.size();
  RemoveOutliers_AngleError(sfm_data_, 2.0);
  const size_t pointcount_angular_filter = sfm_data_.structure.size();
  std::cout << "Outlier removal (remaining #points):\n"
    << "\t initial structure size #3DPoints: " << pointcount_initial << "\n"
    << "\t\t pixel residual filter  #3DPoints: " << pointcount_pixelresidual_filter << "\n"
    << "\t\t angular filter         #3DPoints: " << pointcount_angular_filter << std::endl;

  if (!sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_03_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  // Check that poses & intrinsic cover some measures (after outlier removal)
  const IndexT minPointPerPose = 12; // 6 min
  const IndexT minTrackLength = 3; // 2 min
  if (eraseUnstablePosesAndObservations(sfm_data_, minPointPerPose, minTrackLength))
  {
    // TODO: must ensure that track graph is producing a single connected component

    const size_t pointcount_cleaning = sfm_data_.structure.size();
    std::cout << "Point_cloud cleaning:\n"
      << "\t #3DPoints: " << pointcount_cleaning << "\n";
  }

  // --
  // Final BA. We refine one more time,
  // since some outlier have been removed and so a better solution can be found.
  //--
  const Optimize_Options ba_refine_options(
    ReconstructionEngine::intrinsic_refinement_options_,
    Extrinsic_Parameter_Type::ADJUST_ALL,  // adjust camera motion
    Structure_Parameter_Type::ADJUST_ALL,  // adjust scene structure
    Control_Point_Parameter(),
    this->b_use_motion_prior_);

  b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
  if (b_BA_Status && !sLogging_file_.empty())
  {
    Save(sfm_data_,
      stlplus::create_filespec(stlplus::folder_part(sLogging_file_), "structure_04_outlier_removed", "ply"),
      ESfM_Data(EXTRINSICS | STRUCTURE));
  }

  return b_BA_Status;
}

void GlobalSfMReconstructionEngine_RelativeMotions::Compute_Relative_Rotations
(
  rotation_averaging::RelativeRotations & vec_relatives_R
)
{
  //
  // Build the Relative pose graph from matches:
  //
  /// pairwise view relation between poseIds
  using PoseWiseMatches = std::map< Pair, Pair_Set >;

  // List shared correspondences (pairs) between poses
  PoseWiseMatches poseWiseMatches;
  for (const auto & iterMatches : matches_provider_->pairWise_matches_)
  {
    const Pair pair = iterMatches.first;
    const View * v1 = sfm_data_.GetViews().at(pair.first).get();
    const View * v2 = sfm_data_.GetViews().at(pair.second).get();
    poseWiseMatches[Pair(v1->id_pose, v2->id_pose)].insert(pair);
  }

  C_Progress_display my_progress_bar( poseWiseMatches.size(),
      std::cout, "\n- Relative pose computation -\n" );

#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  // Compute the relative pose from pairwise point matches:
  for (int i = 0; i < static_cast<int>(poseWiseMatches.size()); ++i)
  {
    ++my_progress_bar;
    {
      PoseWiseMatches::const_iterator iter (poseWiseMatches.begin());
      std::advance(iter, i);
      const auto & relative_pose_iterator(*iter);
      const Pair relative_pose_pair = relative_pose_iterator.first;
      const Pair_Set & match_pairs = relative_pose_iterator.second;

      // If a pair has the same ID, discard it
      if (relative_pose_pair.first == relative_pose_pair.second)
      {
        continue;
      }

      // Select common bearing vectors
      if (match_pairs.size() > 1)
      {
        std::cerr << "Compute relative pose between more than two view is not supported" << std::endl;
        continue;
      }

      const Pair pairIterator = *(match_pairs.begin());

      const IndexT I = pairIterator.first;
      const IndexT J = pairIterator.second;

      const View * view_I = sfm_data_.views[I].get();
      const View * view_J = sfm_data_.views[J].get();

      // Check that valid cameras are existing for the pair of view
      if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic) == 0 ||
        sfm_data_.GetIntrinsics().count(view_J->id_intrinsic) == 0)
        continue;


      const IntrinsicBase * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
      const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();

      // Setup corresponding bearing vector
      const matching::IndMatches & matches = matches_provider_->pairWise_matches_.at(pairIterator);
      size_t nbBearing = matches.size();
      Mat x1(2, nbBearing), x2(2, nbBearing);
      nbBearing = 0;
      for (const auto & match : matches)
      {
        x1.col(nbBearing) = ((*cam_I)(cam_I->get_ud_pixel(features_provider_->feats_per_view[I][match.i_].coords().cast<double>()))).hnormalized();
        x2.col(nbBearing++) = ((*cam_J)(cam_J->get_ud_pixel(features_provider_->feats_per_view[J][match.j_].coords().cast<double>()))).hnormalized();
      }

      RelativePose_Info relativePose_info;
      // Compute max authorized error as geometric mean of camera plane tolerated residual error
      relativePose_info.initial_residual_tolerance = std::pow(
        cam_I->imagePlane_toCameraPlaneError(2.5) *
        cam_J->imagePlane_toCameraPlaneError(2.5),
        1./2.);

      // Since we use normalized features, we will use unit image size and intrinsic matrix:
      const std::pair<size_t, size_t> imageSize(1., 1.);
      const Mat3 K  = Mat3::Identity();

      if (!robustRelativePose(K, K, x1, x2, relativePose_info, imageSize, imageSize, 256))
      {
        continue;
      }
      const bool bRefine_using_BA = true;
      if (bRefine_using_BA)
      {
        // Refine the defined scene
        SfM_Data tiny_scene;
        tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
        tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
        tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
        tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

        // Init poses
        const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
        const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

        // Init structure
        const Mat34 P1 = cam_I->get_projective_equivalent(Pose_I);
        const Mat34 P2 = cam_J->get_projective_equivalent(Pose_J);
        Landmarks & landmarks = tiny_scene.structure;
        for (Mat::Index k = 0; k < x1.cols(); ++k)
        {
          const Vec2 x1_ = features_provider_->feats_per_view[I][matches[k].i_].coords().cast<double>();
          const Vec2 x2_ = features_provider_->feats_per_view[J][matches[k].j_].coords().cast<double>();
          Vec3 X;
          TriangulateDLT(P1, x1_.homogeneous(), P2, x2_.homogeneous(), &X);
          Observations obs;
          obs[view_I->id_view] = Observation(x1_, matches[k].i_);
          obs[view_J->id_view] = Observation(x2_, matches[k].j_);
          landmarks[k].obs = obs;
          landmarks[k].X = X;
        }
        // - refine only Structure and Rotations & translations (keep intrinsic constant)
        Bundle_Adjustment_Ceres::BA_Ceres_options options(false, false);
        options.linear_solver_type_ = ceres::DENSE_SCHUR;
        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        const Optimize_Options ba_refine_options
          (Intrinsic_Parameter_Type::NONE, // -> Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL);// adjust scene structure
        if (bundle_adjustment_obj.Adjust(tiny_scene, ba_refine_options))
        {
          // --> to debug: save relative pair geometry on disk
          // std::ostringstream os;
          // os << relative_pose_pair.first << "_" << relative_pose_pair.second << ".ply";
          // Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
          //
          const Mat3 R1 = tiny_scene.poses[view_I->id_pose].rotation();
          const Mat3 R2 = tiny_scene.poses[view_J->id_pose].rotation();
          const Vec3 t1 = tiny_scene.poses[view_I->id_pose].translation();
          const Vec3 t2 = tiny_scene.poses[view_J->id_pose].translation();
          // Compute relative motion and save it
          Mat3 Rrel;
          Vec3 trel;
          RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
          // Update found relative pose
          relativePose_info.relativePose = Pose3(Rrel, -Rrel.transpose() * trel);
        }
      }
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
      {
        // Add the relative rotation to the relative 'rotation' pose graph
        using namespace openMVG::rotation_averaging;
          vec_relatives_R.emplace_back(
            relative_pose_pair.first, relative_pose_pair.second,
            relativePose_info.relativePose.rotation(),
            1.f);
      }
    }
  } // for all relative pose

  // Log input graph to the HTML report
  if (!sLogging_file_.empty() && !sOut_directory_.empty())
  {
    // Log a relative view graph
    {
      std::set<IndexT> set_ViewIds;
      std::transform(sfm_data_.GetViews().begin(), sfm_data_.GetViews().end(),
        std::inserter(set_ViewIds, set_ViewIds.begin()), stl::RetrieveKey());
      graph::indexedGraph putativeGraph(set_ViewIds, getPairs(matches_provider_->pairWise_matches_));
      graph::exportToGraphvizData(
        stlplus::create_filespec(sOut_directory_, "global_relative_rotation_view_graph"),
        putativeGraph);
    }

    // Log a relative pose graph
    {
      std::set<IndexT> set_pose_ids;
      Pair_Set relative_pose_pairs;
      for (const auto & relative_R : vec_relatives_R)
      {
        const Pair relative_pose_indices(relative_R.i, relative_R.j);
        relative_pose_pairs.insert(relative_pose_indices);
        set_pose_ids.insert(relative_R.i);
        set_pose_ids.insert(relative_R.j);
      }
      const std::string sGraph_name = "global_relative_rotation_pose_graph";
      graph::indexedGraph putativeGraph(set_pose_ids, relative_pose_pairs);
      graph::exportToGraphvizData(
        stlplus::create_filespec(sOut_directory_, sGraph_name),
        putativeGraph);
      using namespace htmlDocument;
      std::ostringstream os;

      os << "<br>" << "global_relative_rotation_pose_graph" << "<br>"
         << "<img src=\""
         << stlplus::create_filespec(sOut_directory_, "global_relative_rotation_pose_graph", "svg")
         << "\" height=\"600\">\n";

      html_doc_stream_->pushInfo(os.str());
    }
  }
}

} // namespace sfm
} // namespace openMVG
