/*
* Copyright 2025 Noblis, Inc.
* https://fingerprint.nist.gov/openlqm
*
* Licensed under the Apache License, Version 2.0 (the "License"); you may not
* use this file except in compliance with the License. You may obtain a copy of
* the License at https://www.apache.org/licenses/LICENSE-2.0.
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
* License for the specific language governing permissions and limitations under
* the License.
*
* OpenLQM was developed by Noblis, Inc. under contract to the National
* Institute of Standards and Technology in 2024-2025, as a port of "LQMetric,"
* which was developed by Noblis, Inc. under contract to the Federal Bureau of
* Investigation's Criminal Justice Information Services Division in 2012-2014.
*/

#include <openlqm/openlqm_clarity.hpp>
#include <openlqm/openlqm_impl.hpp>
#include <opencv2/opencv.hpp>

#include <cmath>

namespace OpenLQM {
	namespace Core {
		void GenerateNewClarityMapRandomForest(
			int mapWidth, int mapHeight,
			cv::Mat& clarityMap,
			const unsigned char* pGrayRangeMap,
			const unsigned char* pCurvatureMap,
			const unsigned char* pDirChangeMap,
			const double* pLowFreqMagMap,
			const double* pMaxMagMap,
			const unsigned char* pGrayMedianMap,
			const unsigned char* pGrayCountMap,
			const unsigned char* pValidNeighbors,
			const double* pNormMagMap,
			const unsigned char* pQualMap
		) {
			cv::AutoBuffer<int> votes(static_cast<std::size_t>(OpenLQM::Impl::numClasses));

			cv::Mat m_curvatureMap(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(reinterpret_cast<const void*>(pCurvatureMap)), static_cast<std::size_t>(mapWidth));
			cv::Mat m_curvatureMapA(mapHeight, mapWidth, CV_32FC1);
			m_curvatureMapA.setTo(0);
			absDeviation(m_curvatureMap, m_curvatureMapA, 7);

			cv::Mat m_dirChangeMap(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(reinterpret_cast<const void*>(pDirChangeMap)), static_cast<std::size_t>(mapWidth));
			cv::Mat m_dirChangeMapA(mapHeight, mapWidth, CV_32FC1);
			m_dirChangeMapA.setTo(0);
			absDeviation(m_dirChangeMap, m_dirChangeMapA, 7);

			cv::Mat m_grayCountMap(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(reinterpret_cast<const void*>(pGrayCountMap)), static_cast<std::size_t>(mapWidth));
			cv::Mat m_grayCountMapA(mapHeight, mapWidth, CV_32FC1);
			m_grayCountMapA.setTo(0);
			absDeviation(m_grayCountMap, m_grayCountMapA, 7);

			cv::Mat m_lowFreqMagMap(mapHeight, mapWidth, CV_64FC1, const_cast<void*>(reinterpret_cast<const void*>(pLowFreqMagMap)), static_cast<std::size_t>(mapWidth * 8));
			cv::Mat m_lowFreqMagMapA(mapHeight, mapWidth, CV_32FC1);
			m_lowFreqMagMapA.setTo(0);
			absDeviation(m_lowFreqMagMap, m_lowFreqMagMapA, 7);

			cv::Mat m_maxMagMap(mapHeight, mapWidth, CV_64FC1, const_cast<void*>(reinterpret_cast<const void*>(pMaxMagMap)), static_cast<std::size_t>(mapWidth * 8));
			cv::Mat m_maxMagMapA(mapHeight, mapWidth, CV_32FC1);
			m_maxMagMapA.setTo(0);
			absDeviation(m_maxMagMap, m_maxMagMapA, 7);

			cv::Mat m_grayMedianMap(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(reinterpret_cast<const void*>(pGrayMedianMap)), static_cast<std::size_t>(mapWidth));
			cv::Mat m_grayMedianMapA(mapHeight, mapWidth, CV_32FC1);
			m_grayMedianMapA.setTo(0);
			absDeviation(m_grayMedianMap, m_grayMedianMapA, 7);

			cv::Mat m_validNeighborsMap(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(reinterpret_cast<const void*>(pValidNeighbors)), static_cast<std::size_t>(mapWidth));
			cv::Mat m_validNeighborsMapA(mapHeight, mapWidth, CV_32FC1);
			m_validNeighborsMapA.setTo(0);
			absDeviation(m_validNeighborsMap, m_validNeighborsMapA, 7);

			cv::Mat m_normMagMap(mapHeight, mapWidth, CV_64FC1, const_cast<void*>(reinterpret_cast<const void*>(pNormMagMap)), static_cast<std::size_t>(mapWidth * 8));
			cv::Mat m_normMagMapA(mapHeight, mapWidth, CV_32FC1);
			m_normMagMapA.setTo(0);
			absDeviation(m_normMagMap, m_normMagMapA, 7);

			cv::Mat m_nistQualMap(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(reinterpret_cast<const void*>(pQualMap)), static_cast<std::size_t>(mapWidth));
			cv::Mat m_nistQualMapA(mapHeight, mapWidth, CV_32FC1);
			m_nistQualMapA.setTo(0);
			absDeviation(m_nistQualMap, m_nistQualMapA, 7);

			cv::Mat m_grayRangeMap(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(reinterpret_cast<const void*>(pGrayRangeMap)), static_cast<std::size_t>(mapWidth));
			cv::Mat m_grayRangeMapA(mapHeight, mapWidth, CV_32FC1);
			m_grayRangeMapA.setTo(0);
			absDeviation(m_grayRangeMap, m_grayRangeMapA, 7);

			cv::Mat test_sample(1, LQM_ATTRIBUTES_PER_SAMPLE, CV_32FC1);

			unsigned char* pClarityMap = reinterpret_cast<unsigned char*>(clarityMap.data);
			for (int i=0; i < clarityMap.rows; i++)
			{
				for (int k = 0; k < clarityMap.cols; k++)
				{
					test_sample.at<float>(0,0) = static_cast<float>(std::min(std::min(std::min(i,k),std::min(clarityMap.cols - k -1, clarityMap.rows - i - 1)),15));

					test_sample.at<float>(0,1) = static_cast<float>(m_grayCountMapA.at<float>(i,k));
					test_sample.at<float>(0,2) = static_cast<float>(m_grayCountMap.at<unsigned char>(i,k));

					test_sample.at<float>(0,3) = static_cast<float>(m_curvatureMapA.at<float>(i,k));
					test_sample.at<float>(0,4) = static_cast<float>(m_curvatureMap.at<unsigned char>(i,k));

					test_sample.at<float>(0,5) = static_cast<float>(m_dirChangeMapA.at<float>(i,k));
					test_sample.at<float>(0,6) = static_cast<float>(m_dirChangeMap.at<unsigned char>(i,k));

					test_sample.at<float>(0,7) = static_cast<float>(m_lowFreqMagMapA.at<float>(i,k));
					test_sample.at<float>(0,8) = static_cast<float>(m_lowFreqMagMap.at<double>(i,k));

					test_sample.at<float>(0,9) = static_cast<float>(m_maxMagMapA.at<float>(i,k));
					test_sample.at<float>(0,10) = static_cast<float>(m_maxMagMap.at<double>(i,k));

					test_sample.at<float>(0,11) = static_cast<float>(m_grayMedianMapA.at<float>(i,k));
					test_sample.at<float>(0,12) = static_cast<float>(m_grayMedianMap.at<unsigned char>(i,k));

					test_sample.at<float>(0,13) = static_cast<float>(m_validNeighborsMapA.at<float>(i,k));
					test_sample.at<float>(0,14) = static_cast<float>(m_validNeighborsMap.at<unsigned char>(i,k));

					test_sample.at<float>(0,15) = static_cast<float>(m_normMagMapA.at<float>(i,k));
					test_sample.at<float>(0,16) = static_cast<float>(m_normMagMap.at<double>(i,k));

					test_sample.at<float>(0,17) = static_cast<float>(m_nistQualMapA.at<float>(i,k));
					test_sample.at<float>(0,18) = static_cast<float>(m_nistQualMap.at<unsigned char>(i,k));

					test_sample.at<float>(0,19) = static_cast<float>(m_grayRangeMapA.at<float>(i,k));
					test_sample.at<float>(0,20) = static_cast<float>(m_grayRangeMap.at<unsigned char>(i,k));

					float result = RandomForestPredictClass(OpenLQM::Impl::roots, OpenLQM::Impl::nodes, OpenLQM::Impl::splits, OpenLQM::Impl::numTrees, OpenLQM::Impl::numClasses, votes, test_sample);
					*pClarityMap++ = static_cast<unsigned char>(result);
					// clarityMap.at<unsigned char>(i,k) = (int)result;
				}
			}
		}

		void absDeviation(const cv::Mat& inputFeatureMap, cv::Mat& outputAggFeatureMap, const int filterSize)
		{
			float reducedFilterSize = static_cast<float>(filterSize/2);
			int iNumPadding = static_cast<int>(std::floor(reducedFilterSize));
			int mapWidth = inputFeatureMap.cols;
			int mapHeight = inputFeatureMap.rows;
			int m_padWidth = mapWidth + (iNumPadding * 2);
			int m_padHeight = mapHeight + (iNumPadding * 2);

			cv::Mat m_paddingData(m_padHeight,m_padWidth,inputFeatureMap.type());
			m_paddingData.setTo(0);

			cv::Mat sampleData(filterSize,filterSize,inputFeatureMap.type());
			sampleData.setTo(0);

			int idx2;
			int idx1;

			if (inputFeatureMap.type() == CV_8UC1)
			{
				for (int i = 0; i < mapHeight; ++i)
				{
					for (int k =0; k < mapWidth; ++k)
					{
						idx1 = ((i + iNumPadding) * m_padWidth) + iNumPadding + k;
						idx2 = (i * mapWidth) + k;


						m_paddingData.data[idx1]= inputFeatureMap.data[idx2];
					}

				}

				for (int i = 0; i < iNumPadding; ++i)
				{
					for (int k = 0; k < mapWidth; ++k)
					{
						idx1 = (i * m_padWidth) + iNumPadding + k;
						idx2 = k;
						m_paddingData.data[idx1] = inputFeatureMap.data[idx2];

						idx1 = ((i+iNumPadding + mapHeight) * m_padWidth) + iNumPadding + k;
						idx2 = ((mapHeight - 1) * mapWidth) + k;
						m_paddingData.data[idx1] = inputFeatureMap.data[idx2];
					}
				}


				for (int i = 0; i < m_padHeight; ++i)
				{
					for (int k = 0; k < iNumPadding; ++k)
					{
						idx1 = (i * m_padWidth) + k;
						idx2 = (i * m_padWidth) + iNumPadding;
						m_paddingData.data[idx1] = m_paddingData.data[idx2];

						idx1 = (i * m_padWidth) + iNumPadding + mapWidth + k;
						idx2 = (i * m_padWidth) + iNumPadding + mapWidth - 1;
						m_paddingData.data[idx1] = m_paddingData.data[idx2];
					}
				}
			}
			else
			{
				cv::Mat roiImg = m_paddingData(cv::Rect(iNumPadding,iNumPadding,inputFeatureMap.cols,inputFeatureMap.rows));
				inputFeatureMap.copyTo(roiImg);

				if (inputFeatureMap.type() == CV_32FC1)
				{
					for(int p=0;p<3; p++)
					{
						for(int f=3;f<m_paddingData.cols-3;f++)
						{
							m_paddingData.at<float>(p,f) = inputFeatureMap.at<float>(0,f-3);
						}
					}

					for(int p=m_paddingData.rows -3;p<m_paddingData.rows; p++)
					{
						for(int f=3;f<m_paddingData.cols-3;f++)
						{
							m_paddingData.at<float>(p,f) = inputFeatureMap.at<float>(inputFeatureMap.rows-1,f-3);
						}
					}

					for(int p=0;p<m_paddingData.rows; p++)
					{
						for(int f=0;f<3;f++)
						{
							m_paddingData.at<float>(p,f) = m_paddingData.at<float>(p,3);
						}
					}

					for(int p=0;p<m_paddingData.rows; p++)
					{
						for(int f=m_paddingData.cols-3;f<m_paddingData.cols;f++)
						{
							m_paddingData.at<float>(p,f) = m_paddingData.at<float>(p,m_paddingData.cols-4);
						}
					}


				}
				else
				{
					for(int p=0;p<3; p++)
					{
						for(int f=3;f<m_paddingData.cols-3;f++)
						{
							m_paddingData.at<double>(p,f) = inputFeatureMap.at<double>(p,f-3);
						}
					}
					for(int p=m_paddingData.rows -3;p<m_paddingData.rows; p++)
					{
						for(int f=3;f<m_paddingData.cols-3;f++)
						{
							m_paddingData.at<double>(p,f) = inputFeatureMap.at<double>(inputFeatureMap.rows-1,f-3);
						}
					}

					for(int p=0;p<m_paddingData.rows; p++)
					{
						for(int f=0;f<3;f++)
						{
							m_paddingData.at<double>(p,f) = m_paddingData.at<double>(p,3);
						}
					}

					for(int p=0;p<m_paddingData.rows; p++)
					{
						for(int f=m_paddingData.cols-3;f<m_paddingData.cols;f++)
						{
							m_paddingData.at<double>(p,f) = m_paddingData.at<double>(p,m_paddingData.cols-4);
						}
					}
				}
			}

			for (int y=0; y < m_padHeight - filterSize+1; y++)
			{
				for (int x = 0; x < m_padWidth - filterSize+1; x++)
				{
					for (int i = 0; i < filterSize; i++)
					{
						for (int k = 0; k < filterSize; k++)
						{
							if (inputFeatureMap.type() == CV_8UC1)
							{
								sampleData.at<unsigned char>(i,k) = m_paddingData.at<unsigned char>(y+i,x+k);
							}
							else if (inputFeatureMap.type() == CV_32FC1)
							{
								sampleData.at<float>(i,k) = m_paddingData.at<float>(y+i,x+k);
							}
							else
							{
								sampleData.at<double>(i,k) = m_paddingData.at<double>(y+i,x+k);
							}
						}
					}

					outputAggFeatureMap.at<float>(y,x) = static_cast<float>(calcAbsDeviation(sampleData, filterSize));
				}
			}
		}

		double calcAbsDeviation(const cv::Mat& sData, const int filterSize) {
			double nMidVal;

			if (sData.type() == CV_8UC1)
			{
				nMidVal = static_cast<double>(sData.at<unsigned char>(3,3));
			}
			else if (sData.type() == CV_32FC1)
			{
				nMidVal = static_cast<double>(sData.at<float>(3,3));
			}
			else
			{
				nMidVal = sData.at<double>(3,3);
			}

			double nSum = 0.0;
			double nResult = 0.0;
			double nDiff = 0.0;

			for (int i=0; i < filterSize; i++)
			{
				for (int k = 0; k < filterSize; k++)
				{
					if (i != 3 || k != 3)
					{
						if (sData.type() == CV_8UC1)
						{
							nDiff = nMidVal - static_cast<double>(sData.at<unsigned char>(i,k));
						}
						else if (sData.type() == CV_32FC1)
						{
							nDiff = nMidVal - static_cast<double>(sData.at<float>(i,k));
						}
						else
						{
							nDiff = nMidVal - sData.at<double>(i,k);
						}



						nDiff = std::abs(nDiff);
						nSum +=nDiff;
					}
				}
			}

			nResult = nSum/((filterSize * filterSize)-1);

			return nResult;
		}

		void CropToLargestContiguousRegionOfRidgeFlow(cv::Mat& clarityMap) {
			cv::Mat maskMap;

			clarityMap.copyTo(maskMap);
			int regionCount = 0;
			for (int i = 0; i < maskMap.rows; i++) {
				for (int k = 0; k < maskMap.cols; k++) {
					if (maskMap.at<unsigned char>(i,k) >= 2) {
						maskMap.at<unsigned char>(i,k) = 255;
						regionCount++;
					}
				}
			}

			if (regionCount < 100) {
				regionCount = 0;
				for (int i = 0; i < maskMap.rows; i++) {
					for (int k = 0; k < maskMap.cols; k++) {
						if (maskMap.at<unsigned char>(i,k) >= 1) {
							maskMap.at<unsigned char>(i,k) = 255;
							regionCount++;
						} else {
							maskMap.at<unsigned char>(i,k) = 0;
						}

					}
				}
			} else {
				for (int i = 0; i < maskMap.rows; i++) {
					for (int k = 0; k < maskMap.cols; k++) {
						if (maskMap.at<unsigned char>(i,k) < 2) {
							maskMap.at<unsigned char>(i,k) = 0;
						}
					}
				}
			}

			if (regionCount < 100) {
				clarityMap.setTo(0);
				return;
			}

			cv::Mat smoothedMap;
			cv::blur(clarityMap,smoothedMap,cv::Size(5,5));

			cv::Mat const shape = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3, 3));
			cv::erode(maskMap,maskMap,shape);

			cv::blur(maskMap,maskMap,cv::Size(3,3));

			FillContourHoles(maskMap);

			std::vector<std::vector<cv::Point> > contourPoints;
			std::vector<cv::Vec4i> hierarchy;

			cv::findContours( maskMap, contourPoints, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));
			if (contourPoints.empty())
			        return;

			//Find the largest contour
			int maxContourSizeIndex = 0;
			int maxSize = 0;
			for (int i = 0; i < static_cast<int>(contourPoints.size()); i++) {
				if (static_cast<int>(contourPoints[static_cast<std::vector<std::vector<cv::Point>>::size_type>(i)].size()) > maxSize) {
					maxSize = static_cast<int>(contourPoints[static_cast<std::vector<std::vector<cv::Point>>::size_type>(i)].size());
					maxContourSizeIndex = i;
				}
			}

			maskMap.setTo(0);
			//Check Boundary
			for (int i = 0; i < static_cast<int>(contourPoints[static_cast<std::size_t>(maxContourSizeIndex)].size()); i++) {
				maskMap.at<unsigned char>(contourPoints[static_cast<std::vector<std::vector<cv::Point>>::size_type>(maxContourSizeIndex)][static_cast<std::vector<cv::Point>::size_type>(i)].y, contourPoints[static_cast<std::vector<std::vector<cv::Point>>::size_type>(maxContourSizeIndex)][static_cast<std::vector<cv::Point>::size_type>(i)].x) = 255;
			}

			std::vector<std::vector<cv::Point>> contourToBeFilled;
			contourToBeFilled.push_back(contourPoints[static_cast<std::vector<std::vector<cv::Point>>::size_type>(maxContourSizeIndex)]);
			cv::fillPoly(maskMap,contourToBeFilled, cv::Scalar(1),8);

			cv::multiply(maskMap, smoothedMap, clarityMap);
		}

		void CropToRoi(cv::Mat& clarityMap, const std::vector<OpenLQM::Coordinate>& roi) {
			cv::Mat maskMap = cv::Mat::zeros(clarityMap.rows, clarityMap.cols, clarityMap.type());

			cv::Mat smoothedMap;
			cv::blur(clarityMap,smoothedMap, cv::Size(5,5));

			std::vector<std::vector<cv::Point> > contourPoints;
			contourPoints.resize(1);
			std::vector<cv::Point>& contour = contourPoints[0];
			for (const OpenLQM::Coordinate& coord : roi) {
				contour.emplace_back(coord.x, coord.y);
			}

			cv::fillPoly(maskMap, contourPoints, cv::Scalar(1),8);

			cv::multiply(maskMap, smoothedMap, clarityMap);
		}

		void FillContourHoles(cv::Mat& input) {
			cv::Mat holes = input.clone();
			cv::floodFill(holes, cv::Point2i(0, 0), cv::Scalar(1));
			for (int i = 0; i < input.rows * input.cols; ++i) {
				if (holes.data[i] == 0) {
					input.data[i] = 255;
				}
			}
		}

		float RandomForestPredictClass(const std::vector<int>& roots, const std::vector<cv::ml::DTrees::Node>& nodes, const std::vector<cv::ml::DTrees::Split>& splits, const int numTrees, const int numClasses, cv::AutoBuffer<int>& votes, const cv::Mat& sample) {
			const float* pSample = reinterpret_cast<const float*>(sample.data);

			double result = -1;
			int max_nvotes = 0;
			int* pVotes = votes;
			memset(pVotes, 0, sizeof(*pVotes)*static_cast<std::size_t>(numClasses));

			for (int k = 0; k < numTrees; ++k) {
				// Get the predicted node
				int root = roots[static_cast<std::size_t>(k)];
				const cv::ml::DTrees::Node* node = &nodes[static_cast<std::size_t>(root)];

				while (node->left != -1) {
					int splitIndex = node->split;
					int dir = 0;
					const cv::ml::DTrees::Split* pSplit = nullptr;
					for (; !dir && splitIndex >= 0; splitIndex = pSplit->next) {
						pSplit = &splits[static_cast<std::size_t>(splitIndex)];
						float val = pSample[static_cast<std::size_t>(pSplit->varIdx)];
						dir = (val <= pSplit->c ? -1 : 1);

						if (pSplit->inversed) {
							dir = -dir;
						}
					}

					if (!dir) {
						dir = node->defaultDir;
					}

					int nextNodeIndex = (dir < 0 ? node->left : node->right);
					node = &nodes[static_cast<std::size_t>(nextNodeIndex)];
				}

				// Got the predicted node. Update the vote record and current winner.
				int nvotes = ++pVotes[node->classIdx];
				if (nvotes > max_nvotes) {
					max_nvotes = nvotes;
					result = node->value;
				}
			}

			return static_cast<float>(result);
		}
	}
}
