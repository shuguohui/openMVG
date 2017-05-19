// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
#define OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/sfm_engine.hpp"

namespace htmlDocument { class htmlDocumentStream; }
namespace openMVG { namespace matching { struct PairWiseMatches; } }
namespace openMVG { namespace sfm { struct Features_Provider; } }
namespace openMVG { namespace sfm { struct Matches_Provider; } }
namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG{
namespace sfm{

/// Global SfM Pipeline Reconstruction Engine.
/// - Method: Global Fusion of Relative Motions.
class GlobalSfMReconstructionEngine_RelativeMotions : public ReconstructionEngine
{
public:
	// Track statistics are the track length and mean reprojection error.
	typedef std::pair<int, double> TrackStatistics;
	typedef std::pair<IndexT, TrackStatistics> GridCellElement;
	typedef std::unordered_map<int64_t, std::vector<GridCellElement> > ImageGrid;

public:

  GlobalSfMReconstructionEngine_RelativeMotions(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_RelativeMotions() override;

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);

  void SetRotationAveragingMethod(ERotationAveragingMethod eRotationAveragingMethod);
  void SetTranslationAveragingMethod(ETranslationAveragingMethod eTranslation_averaging_method_);

  bool Process() override;

protected:
  /// Compute from relative rotations the global rotations of the camera poses
  bool Compute_Global_Rotations
  (
    const openMVG::rotation_averaging::RelativeRotations & vec_relatives_R,
    Hash_Map<IndexT, Mat3> & map_globalR
  );

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations
  (
    const Hash_Map<IndexT, Mat3> & global_rotations,
    matching::PairWiseMatches & tripletWise_matches
  );

  /// Compute the initial structure of the scene
  bool Compute_Initial_Structure
  (
    matching::PairWiseMatches & tripletWise_matches
  );

  bool SelectGoodTracksForBundleAdjustment();

  void ComputeTrackStatistics(const int long_track_length_threshold,
							 std::unordered_map<IndexT, TrackStatistics>& track_statistics);

  void SelectBestTracksFromEachImageGridCell(IndexT viewID, int image_grid_cell_size,
											const std::unordered_set<IndexT>& trackids,
											const std::unordered_map<IndexT, TrackStatistics>& track_statistics, 
											std::unordered_set<IndexT>& tracks_to_optimize);

  void SelectTopRankedTracksInView(const std::unordered_map<IndexT, TrackStatistics>& track_statistics, 
									const std::unordered_set<IndexT>& trackids,
									IndexT viewID, int min_num_optimized_tracks_per_view, 
									std::unordered_set<IndexT>& tracks_to_optimize);
  // Adjust the scene (& remove outliers)
  bool Adjust();

private:
  /// Compute relative rotations
  void Compute_Relative_Rotations
  (
    openMVG::rotation_averaging::RelativeRotations & vec_relatives_R
  );

  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> html_doc_stream_;
  std::string sLogging_file_;

  // Parameter
  ERotationAveragingMethod eRotation_averaging_method_;
  ETranslationAveragingMethod eTranslation_averaging_method_;

  //-- Data provider
  Features_Provider  * features_provider_;
  Matches_Provider  * matches_provider_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
