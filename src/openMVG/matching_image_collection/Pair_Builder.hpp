// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IMAGE_COLLECTION_PAIR_BUILDER_HPP
#define OPENMVG_MATCHING_IMAGE_COLLECTION_PAIR_BUILDER_HPP

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "openMVG/types.hpp"
#include "openMVG/stl/split.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/features/regions.hpp"
#include "third_party/libvot/vocab_tree.h"
#include "third_party/libvot/vot_pipeline.h"
#include "third_party/libvot/io_utils.h"

#define MAX_ARRAY_SIZE 8388608  // 2^23
namespace openMVG {

	inline std::vector<size_t> RandomSample(size_t total, size_t sample_num)
	{
		std::vector<size_t> total_numbers(total, 0);
		std::iota(total_numbers.begin(), total_numbers.end(), 0);
		std::random_shuffle(total_numbers.begin(), total_numbers.end());
		std::vector<size_t> samples(total_numbers.begin(), total_numbers.begin() + sample_num);

		return samples;
	}

	inline Pair_Set getPairsFromVocabTree(const sfm::SfM_Data& sfm_data, const std::shared_ptr<sfm::Regions_Provider>& regions_provider)
	{
		Pair_Set pairs;
		int depth = 6;
		int branch_num = 8;
		size_t num_view = sfm_data.views.size();
		size_t memory_size = tw::IO::GetAvailMem() / (1024 * 1024);	// convert to mb
		size_t tree_memory_size = FDIM * sizeof(DTYPE) * (size_t)pow((double)branch_num, (double)(depth + 1)) / (1024 * 1024);
		size_t max_siftfile_num = (memory_size - tree_memory_size) / 2;
		size_t sample_size = num_view > max_siftfile_num ? max_siftfile_num : num_view;
		std::vector<size_t> siftfile_samples = RandomSample(num_view, sample_size);

		size_t total_keys = 0;
		std::vector<vot::SiftData> sift_data;
		sift_data.resize(sample_size);
		const sfm::Regions_Provider* provider = regions_provider.get();
		for (size_t i = 0; i < sample_size; i++)
		{
			total_keys += provider->get(siftfile_samples[i])->RegionCount();
		}

		std::cout << "[Build Tree] Total sift keys  " << total_keys << std::endl;
		if (total_keys == 0)
		{
			std::cerr << "[Build Tree] Error: No sift keys input, maybe the sift type is wrong. Exit...\n";
			return pairs;
		}
		
		const size_t len = (size_t)total_keys * FDIM;
		const size_t num_arrays = len / MAX_ARRAY_SIZE + ((len % MAX_ARRAY_SIZE) == 0 ? 0 : 1);

		// allocate a big chunk of memory to each array
		std::cout << "[Build Tree] Allocate " << len << " bytes memory into " << num_arrays << " arrays\n";
		DTYPE **mem = new DTYPE *[num_arrays];
		size_t remain_length = len;
		for (size_t i = 0; i < num_arrays; i++) 
		{
			const size_t len_curr = remain_length > MAX_ARRAY_SIZE ? MAX_ARRAY_SIZE : remain_length;
			mem[i] = new DTYPE[len_curr];
			remain_length -= len_curr;
		}
		assert(remain_length == 0);
		
		// allocate a pointer array to sift data
		DTYPE **mem_pointer = new DTYPE *[total_keys];
		size_t off = 0;
		size_t curr_key = 0;
		int curr_array = 0;
		for (size_t i = 0; i < sample_size; i++)
		{
			features::Regions* region = provider->get(siftfile_samples[i]).get();
			int num_keys = region->RegionCount();
			if (num_keys > 0) 
			{
				DTYPE *dp = (DTYPE*)region->DescriptorRawData();
				for (int j = 0; j < num_keys; j++) 
				{
					for (int k = 0; k < FDIM; k++)
					{
						mem[curr_array][off + k] = dp[j * FDIM + k];
					}
					mem_pointer[curr_key] = mem[curr_array] + off;
					curr_key++;
					off += FDIM;
					if (off == MAX_ARRAY_SIZE)
					{
						off = 0;
						curr_array++;
					}
				}
			}
		}
		
		// build a vocabulary tree using sift keys
#ifdef OPENMVG_USE_OPENMP
		int thread_num = omp_get_max_threads();
#else 
		int thread_num = 1;
#endif
		vot::VocabTree vt;
		if (!vt.BuildTree(total_keys, FDIM, depth, branch_num, mem_pointer, thread_num))
		{
			return pairs;
		}
		vt.SetConstantWeight();
		for (size_t i = 0; i < num_view; ++i)
		{
			features::Regions* region = provider->get(i).get();

			const double mag = vt.AddImage2Tree( i, region->RegionCount(), (DTYPE*)region->DescriptorRawData(), thread_num);
			std::cout << "[BuildDB] Add image #" <<  i << " to database\n";
		}
		vt.ComputeTFIDFWeight(num_view);
		vt.NormalizeDatabase(0, num_view);

		size_t db_image_num = vt.database_image_num;
		size_t num_matches = 100;
#ifdef OPENMVG_USE_OPENMP
		std::vector<std::vector<float> > scores(thread_num);
		std::vector<std::vector<size_t> > indexed_scores(thread_num);
		for (int i = 0; i < thread_num; i++) 
		{
			scores[i].resize(db_image_num);
			indexed_scores[i].resize(db_image_num);
			std::iota(indexed_scores[i].begin(), indexed_scores[i].end(), 0);
		}
#pragma omp parallel for schedule(dynamic)
#else 
		std::vector<float> scores(db_image_num);
		std::vector<size_t> indexed_scores(db_image_num);
		std::iota(indexed_scores.begin(), indexed_scores.end(), 0);
#endif
		for (int i = 0; i < num_view; ++i)
		{
			features::Regions* region = provider->get(i).get();
#ifdef OPENMVG_USE_OPENMP
			int thread_num = omp_get_thread_num();
			std::vector<float>& vscores = scores[thread_num];
			std::vector<size_t>& vindexed_scores = indexed_scores[thread_num];
			vt.Query(region->RegionCount(), (DTYPE*)region->DescriptorRawData(), &vscores[0]);
			std::sort(vindexed_scores.begin(), vindexed_scores.end(),[&](size_t i0, size_t i1) {return vscores[i0] > vscores[i1]; });
			size_t num_length = std::min(num_matches, vindexed_scores.size());
#pragma omp critical 
			{
				for (size_t j = 0; j < num_length; j++)
				{
					if (i == vindexed_scores[j])
						continue;
					else if (i >vindexed_scores[j])
						pairs.insert(std::make_pair(vindexed_scores[j], i));
					else
						pairs.insert(std::make_pair(i, vindexed_scores[j]));
				}
			}
			
#else
			vt.Query(region->RegionCount(), (DTYPE*)region->DescriptorRawData(), &scores[0]);

			std::sort(indexed_scores.begin(), indexed_scores.end(),
				[&](size_t i0, size_t i1) {return scores[i0] > scores[i1]; });
			size_t num_length = std::min(num_matches, indexed_scores.size());
			for (size_t j = 0; j < num_length; j++)
			{
				if (i == indexed_scores[j])
					continue;
				else if (i >indexed_scores[j])
					pairs.insert(std::make_pair(indexed_scores[j], i));
				else
					pairs.insert(std::make_pair(i, indexed_scores[j]));
			}
#endif
		}
		vt.ClearTree();

		return pairs;
	}
	/// Generate all the (I,J) pairs of the upper diagonal of the NxN matrix
	inline Pair_Set exhaustivePairs(const size_t N)
	{
		Pair_Set pairs;
		for (IndexT I = 0; I < static_cast<IndexT>(N); ++I)
			for (IndexT J = I + 1; J < static_cast<IndexT>(N); ++J)
				pairs.insert(std::make_pair(I, J));

		return pairs;
	}

	/// Generate the pairs that have a distance inferior to the overlapSize
	/// Usable to match video sequence
	inline Pair_Set contiguousWithOverlap(const size_t N, const size_t overlapSize)
	{
		Pair_Set pairs;
		for (IndexT I = 0; I < static_cast<IndexT>(N); ++I)
			for (IndexT J = I + 1; J < I + 1 + overlapSize && J < static_cast<IndexT>(N); ++J)
				pairs.insert(std::make_pair(I, J));
		return pairs;
	}

	/// Load a set of Pair_Set from a file
	/// I J K L (pair that link I)
	inline bool loadPairs(
		const size_t N,  // number of image in the current project (to check index validity)
		const std::string &sFileName, // filename of the list file,
		Pair_Set & pairs)  // output pairs read from the list file
	{
		std::ifstream in(sFileName.c_str());
		if (!in.is_open())
		{
			std::cerr << std::endl
				<< "loadPairs: Impossible to read the specified file: \"" << sFileName << "\"." << std::endl;
			return false;
		}
		std::string sValue;
		std::vector<std::string> vec_str;
		while (std::getline(in, sValue))
		{
			vec_str.clear();
			stl::split(sValue, ' ', vec_str);
			const IndexT str_size(vec_str.size());
			if (str_size < 2)
			{
				std::cerr << "loadPairs: Invalid input file: \"" << sFileName << "\"." << std::endl;
				return false;
			}
			std::stringstream oss;
			oss.clear(); oss.str(vec_str[0]);
			IndexT I, J;
			oss >> I;
			for (IndexT i = 1; i < str_size; ++i)
			{
				oss.clear(); oss.str(vec_str[i]);
				oss >> J;
				if (I > N - 1 || J > N - 1) //I&J always > 0 since we use unsigned type
				{
					std::cerr << "loadPairs: Invalid input file. Image out of range. "
						<< "I: " << I << " J:" << J << " N:" << N << std::endl
						<< "File: \"" << sFileName << "\"." << std::endl;
					return false;
				}
				if (I == J)
				{
					std::cerr << "loadPairs: Invalid input file. Image " << I << " see itself. File: \"" << sFileName << "\"." << std::endl;
					return false;
				}
				pairs.insert((I < J) ? std::make_pair(I, J) : std::make_pair(J, I));
			}
		}
		in.close();
		return true;
	}

	/// Save a set of Pair_Set to a file (one pair per line)
	/// I J
	/// I K
	/// ...
	inline bool savePairs(const std::string &sFileName, const Pair_Set & pairs)
	{
		std::ofstream outStream(sFileName.c_str());
		if (!outStream.is_open()) {
			std::cerr << std::endl
				<< "savePairs: Impossible to open the output specified file: \"" << sFileName << "\"." << std::endl;
			return false;
		}
		for (const auto & cur_pair : pairs)
		{
			outStream << cur_pair.first << ' ' << cur_pair.second << '\n';
		}
		bool bOk = !outStream.bad();
		outStream.close();
		return bOk;
	}

} // namespace openMVG

#endif // OPENMVG_MATCHING_IMAGE_COLLECTION_PAIR_BUILDER_HPP
