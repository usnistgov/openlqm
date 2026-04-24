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

#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <openlqm/openlqm.hpp>
#include <openlqm/openlqm_minutia.hpp>
#include <openlqm/openlqm_clarity.hpp>

#include <opencv2/opencv.hpp>
#include <opencv2/core/utils/logger.hpp>

#include <stack>
#include <tuple>
#include <stdio.h>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <iostream>
#include <string>
#include <memory>
#include <cmath>
#include <errno.h>

#define ROUND_OUT(x) std::round(x * 10.0) / 10.0
#include <iomanip>
#include <openlqm/openlqm_img_util.hpp>

#include <filesystem>
#ifdef _WIN32
	#include <stdlib.h>
#else
	#include <dlfcn.h>
#endif

namespace OpenLQM {
	namespace Impl {
		std::string GetInstalledModelPath() {
			#ifdef _WIN32
			char* pModel = std::getenv("OPENLQM_PATH");
			if (!pModel) {
				return std::string(OPENLQM_MODEL_INSTALL_PATH);
			}
			std::string olqmPath(pModel);
			return olqmPath + "/openlqm_model.xml";
			#else
			Dl_info info;
			if (dladdr(reinterpret_cast<void*>(&OpenLQM::Supplement::LoadModel), &info)) {
				return std::filesystem::path(info.dli_fname).parent_path().parent_path().string() + "/share/openlqm_model.xml";
			}
			throw std::runtime_error("dladdr failed to find olqm's location");
			#endif
		}

		std::string modelPath;
		cv::Ptr<cv::ml::RTrees> loadedModel = cv::ml::RTrees::load(GetInstalledModelPath());
		std::vector<int> roots = loadedModel->getRoots();
		std::vector<cv::ml::DTrees::Node> nodes = loadedModel->getNodes();
		std::vector<cv::ml::DTrees::Split> splits = loadedModel->getSplits();

		int GetClassCount(const std::vector<cv::ml::DTrees::Node>& nodes_) {
			int maxLabel = -1;
			for (const cv::ml::DTrees::Node& n : nodes_) {
				if (n.left < 0 && n.right < 0) {
					maxLabel = std::max(maxLabel, n.classIdx);
				}
			}

			return maxLabel + 1;
		}

		int numTrees = static_cast<int>(roots.size());
		int numClasses = GetClassCount(nodes);
	}

	void BinarizeMap(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, unsigned char threshold);
	void AreaCellCount(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius);
	void DilateErodeCore(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius, unsigned char val);
	void Dilate(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius = 1);
	void Erode(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius = 1);
	void Open(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius = 1);
	void Close(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius = 1);
	void OpenClose(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius = 1);
	int CountNonZero(const std::vector<unsigned char>& src);
	void MarkRegions(std::vector<unsigned char>& src, int& regionCount, int width, int height); // Input image binarized (0|1) output image (0|2+) output return # of regions
	void RemoveRegions(std::vector<unsigned char>& src, int width, int height, const std::vector<bool>& keepRegions);
	void Mask(std::vector<unsigned char>& target, const std::vector<unsigned char>& mask);
	int SumAllCells(const std::vector<unsigned char>& src);
	int ScaleValues(int inValue, int inMin, int inMax, int outMin, int outMax);

	bool GetAllMetricsAndQualityMapFromFingerprint_Base(const OpenLQM::Fingerprint& inpImg, OpenLQM::Metrics& metrics, OpenLQM::QualityMap& qualMap, bool inputQualityMap, bool justQualityMap, OpenLQM::Supplement::QualityMeasures* pQualMeasures = nullptr, OpenLQM::Supplement::FeatureMaps* pFeatureMaps = nullptr);

	namespace Core {
		constexpr inline int RESULT_VECTOR_LENGTH = 6;
		constexpr inline double PIXEL_MAPPING = 24.22; // 1 SQmm = 24.22 pixels, 1 SQin = 15265 pixels

		struct GlobalInitializer {
			GlobalInitializer() {
				cv::utils::logging::setLogLevel(cv::utils::logging::LOG_LEVEL_SILENT);
			}
		};

		void FloodFill(std::vector<unsigned char>& src, int width, int height, int x, int y, unsigned char targetValue, unsigned char newValue);
		void FloodFillImpl(std::vector<unsigned char>& src, int width, int height, int x, int y, unsigned char targetValue, unsigned char newValue);

		struct ImageAttributes {
			int Width = 0;
			int Height = 0;
			int Resolution = 0;
			int Mean = 0;
			int Quartile[5] = {};
			int Range = 0;
			int GrayValueCount = 0;
			int PixelCount = 0;
			std::array<int, 256> Hist = {};
			bool BasedOnMask = false;

			ImageAttributes();
			ImageAttributes(const std::vector<unsigned char>& buf, int width, int height, int resolution = 0);

			void Init(const std::vector<unsigned char>& buf, int width, int height, int resolution = 0);
			void Load(const std::vector<unsigned char>& buf, int width, int height, int resolution = 0);
			void CalcStatisticsFromHistogram();
			void Recalc(const std::vector<unsigned char>& buf);
		};

		struct AggregateQuality {
			static constexpr int MAP_TO_SQ_IN_FACTOR = 15625; // 15625 (125 * 125) pixels in 1 sq in
			static constexpr double MAP_TO_SQ_MM_FACTOR = 24.22; // 24.22  (4.9213 * 4.9213) pixels in 1 sq mm
			std::array<int, RESULT_VECTOR_LENGTH> m_totalArea = {};
			std::array<int, RESULT_VECTOR_LENGTH> m_totalFilteredArea = {};
			std::array<int, RESULT_VECTOR_LENGTH> m_largestContiguousArea = {};
			std::array<int, RESULT_VECTOR_LENGTH> m_areaWithinGFA = {};
			std::array<int, RESULT_VECTOR_LENGTH> m_areaWithinLCA = {};
			std::array<int, RESULT_VECTOR_LENGTH> m_regionCount = {};
			std::array<int, RESULT_VECTOR_LENGTH> m_consistencyWithinGFA = {};
			std::array<int, RESULT_VECTOR_LENGTH> m_consistencyWithinLCA = {};
			std::array<int, RESULT_VECTOR_LENGTH> m_openCloseRadii = { 0, 0, 2, 1, 1, 1 };
			std::array<int, RESULT_VECTOR_LENGTH> m_ignoreThreshold = { 0, 0, 49, 49, 16, 16 };
			std::map<int, std::string> m_colorNameDictionary{ {0, "black"}, {1, "red"}, {2, "yellow"}, {3, "green"}, {4, "blue"}, {5, "cyan"} };
			std::map<int, std::string> m_colorMeaningDictionary{ {0, "background"}, {1, "impression"}, {2, "ridge flow"}, {3, "good ridge flow"}, {4, "clear level-3 detail"}, {5, "perfect level-3 detail"} };

			std::vector<unsigned char> LargestContiguousAreaOfRidgeFlow = {};
			std::vector<unsigned char> GoodFlowAreaOfRidgeFlow = {};
			int OverallQuality_GFA_Consist = 0;
			int OverallQuality_LCA_Consist = 0;
			int OverallClarity = 0;

			AggregateQuality(const std::vector<unsigned char>& LocQMap, int width, int height, ImageAttributes& ImgAttrib);
			void CalcOverallQuality_Consist();
			int CalcOverallQuality_Consist(std::array<int, RESULT_VECTOR_LENGTH>& ConsistencyArray);
		};

		bool CheckHeuristics(double LCAQ1plus, double LCAQ3plus, int minQ2plus);
		double CalculateLQMetric12(double predictedMatchScore);
		int PredictValueDeterminations(double predictedMatcherScore, bool determFlag);

		long EmbedMinutiaeDetail(
			const LQMParams* pLqmParams,
			unsigned char* grayImage,
			unsigned char* binaryImage,
			int inputResolution,
			int imgWidth, int imgHeight,
			char* DirMap,
			char* QualMap,
			int maxMinutiae,
			MinutiaOut* MinutiaeArray,
			unsigned char* gray10pctMap,  // grayscale of 10th %ile (dark)
			unsigned char* gray90pctMap,  // grayscale of 90th %ile (light)
			unsigned char* grayMedianMap, // median grayscale
			unsigned char* grayRangeMap, // median grayscale
			unsigned char* grayCountMap,   // number of gray values in use
			double* maxMagnitude,
			double* normMagnitude,
			double* lowFreqMagnitude,
			unsigned char* directionChangeMap, // changes in direction map
			unsigned char* validNeighbors,
			unsigned char* curvatureMap
		);

		int ProcessClarityMap(
			unsigned char* pOutClarityMapNonLCA, int mapWidth, int mapHeight,
			const unsigned char* pLocQ,
			const unsigned char* pGrayRangeMap,
			const unsigned char* pCurvatureMap,
			const unsigned char* pDirChangeMap,
			const double* pLowFreqMagMap,
			const double* pMaxMagMap,
			const unsigned char* pGrayMedianMap,
			const unsigned char* pGrayCountMap,
			const unsigned char* pValidNeighbors,
			const double* pNormMagMap,
			const unsigned char* pQualMap,
			const std::vector<OpenLQM::Coordinate>* roi = nullptr
		);

		bool CalculateMetrics(
			const unsigned char* pInputImg, int inputWidth, int inputHeight, int resolution,
			const unsigned char* pLocQ, int mapWidth, int mapHeight,
			MinutiaOut* pMinutiae, int minutiaCount,
			OpenLQM::Metrics& metrics,
			OpenLQM::Supplement::QualityMeasures* pQualMeasures = nullptr
		);

		template<class T>
		void CopyMatToFeatureMap(const cv::Mat& mat, OpenLQM::FeatureMapBase<T>& featureMap) {
			featureMap.width = static_cast<unsigned int>(mat.cols);
			featureMap.height = static_cast<unsigned int>(mat.rows);
			featureMap.buffer.resize(mat.total());
			memcpy(featureMap.buffer.data(), mat.data, mat.total() * sizeof(T));
		}

		template<class T>
		void CopyVectorToFeatureMap(const std::vector<T> vec, unsigned int width, unsigned int height, OpenLQM::FeatureMapBase<T>& featureMap) {
			featureMap.width = width;
			featureMap.height = height;
			featureMap.buffer = vec;
		}
	}
}

OpenLQM::Core::GlobalInitializer globalInit;

namespace {

unsigned int ConvertDpiToDotsPerMeter(double dpi)
{
	return static_cast<unsigned int>(std::lround((dpi * 10000.0) / 254.0));
}

unsigned int ConvertDpcmToDotsPerMeter(double dpcm)
{
	return static_cast<unsigned int>(std::lround(dpcm * 100.0));
}

unsigned int ReadTiffDotsPerMeterX(const std::vector<unsigned char>& data, std::size_t tiffBase, std::size_t tiffEnd)
{
	if (tiffEnd > data.size() || tiffBase + 8 > tiffEnd) {
		return 0;
	}

	bool le = false;
	if (data[tiffBase] == 0x49 && data[tiffBase + 1] == 0x49) {
		le = true;
	} else if (data[tiffBase] == 0x4D && data[tiffBase + 1] == 0x4D) {
		le = false;
	} else {
		return 0;
	}

	auto u16 = [&](std::size_t o) -> uint16_t {
		if (o + 2 > tiffEnd) return 0;
		return le ? static_cast<uint16_t>(data[o] | (static_cast<uint16_t>(data[o+1]) << 8))
		          : static_cast<uint16_t>((static_cast<uint16_t>(data[o]) << 8) | data[o+1]);
	};
	auto u32 = [&](std::size_t o) -> uint32_t {
		if (o + 4 > tiffEnd) return 0;
		return le ? (static_cast<uint32_t>(data[o]) |
		             (static_cast<uint32_t>(data[o+1]) << 8) |
		             (static_cast<uint32_t>(data[o+2]) << 16) |
		             (static_cast<uint32_t>(data[o+3]) << 24))
		          : ((static_cast<uint32_t>(data[o])   << 24) |
		             (static_cast<uint32_t>(data[o+1]) << 16) |
		             (static_cast<uint32_t>(data[o+2]) << 8)  |
		              static_cast<uint32_t>(data[o+3]));
	};

	if (u16(tiffBase + 2) != 42) {
		return 0;
	}

	const uint32_t ifdRelOffset = u32(tiffBase + 4);
	if (ifdRelOffset > tiffEnd - tiffBase) {
		return 0;
	}
	const std::size_t ifdOff = tiffBase + ifdRelOffset;
	if (ifdOff + 2 > tiffEnd) {
		return 0;
	}

	const uint16_t nEntries = u16(ifdOff);
	double xRes = 0.0;
	uint16_t resUnit = 2;
	for (uint16_t i = 0; i < nEntries; ++i) {
		const std::size_t e = ifdOff + 2 + (static_cast<std::size_t>(i) * 12);
		if (e + 12 > tiffEnd) break;
		const uint16_t tag  = u16(e);
		const uint16_t type = u16(e + 2);
		const uint32_t count = u32(e + 4);
		const uint32_t valueOrOffset = u32(e + 8);
		if (tag == 282 && type == 5 && count > 0) { // XResolution RATIONAL
			if (valueOrOffset > tiffEnd - tiffBase) {
				continue;
			}
			const std::size_t roff = tiffBase + valueOrOffset;
			if (roff + 8 > tiffEnd) {
				continue;
			}
			const uint32_t num = u32(roff);
			const uint32_t den = u32(roff + 4);
			if (den) xRes = static_cast<double>(num) / den;
		}
		if (tag == 296 && type == 3 && count > 0) {
			resUnit = u16(e + 8);
		}
	}

	if (xRes > 0.0) {
		if (resUnit == 2) return ConvertDpiToDotsPerMeter(xRes);
		if (resUnit == 3) return ConvertDpcmToDotsPerMeter(xRes);
	}
	return 0;
}

unsigned int ReadDotsPerMeterX(const std::vector<unsigned char>& data)
{
	const std::size_t sz = data.size();
	if (sz < 8) return 0;

	// ---- JPEG (FF D8) ----
	if (data[0] == 0xFF && data[1] == 0xD8) {
		std::size_t pos = 2;
		unsigned int jfifDpm = 0;
		while (pos < sz) {
			if (data[pos] != 0xFF) break;
			while (pos < sz && data[pos] == 0xFF) {
				++pos;
			}
			if (pos >= sz) break;
			const uint8_t marker = data[pos++];
			if (marker == 0xD9 || marker == 0xDA) break; // EOI / SOS
			if (marker == 0x01 || (marker >= 0xD0 && marker <= 0xD7)) continue;
			if (pos + 2 > sz) break;
			const auto segLen = static_cast<uint16_t>(
				(static_cast<uint16_t>(data[pos]) << 8) | data[pos + 1]);
			if (segLen < 2 || pos + segLen > sz) break;
			const std::size_t segmentStart = pos + 2;
			const std::size_t segmentEnd = pos + segLen;

			if (marker == 0xE0 && segmentStart + 14 <= segmentEnd &&
			    data[segmentStart]=='J' && data[segmentStart+1]=='F' && data[segmentStart+2]=='I' &&
			    data[segmentStart+3]=='F' && data[segmentStart+4]==0) {
				const uint8_t units = data[segmentStart + 7];
				const auto xd = static_cast<uint16_t>(
					(static_cast<uint16_t>(data[segmentStart + 8]) << 8) | data[segmentStart + 9]);
				if (xd > 0) {
					if (units == 1) jfifDpm = ConvertDpiToDotsPerMeter(xd);
					if (units == 2) jfifDpm = ConvertDpcmToDotsPerMeter(xd);
				}
			}
			if (marker == 0xE1 && segmentStart + 14 <= segmentEnd &&
			    data[segmentStart]=='E' && data[segmentStart+1]=='x' && data[segmentStart+2]=='i' &&
			    data[segmentStart+3]=='f' && data[segmentStart+4]==0 && data[segmentStart+5]==0) {
				const unsigned int exifDpm = ReadTiffDotsPerMeterX(data, segmentStart + 6, segmentEnd);
				if (exifDpm > 0) return exifDpm;
			}
			pos += segLen;
		}
		return jfifDpm;
	}

	// ---- PNG (89 50 4E 47) ----
	if (data[0]==0x89 && data[1]==0x50 && data[2]==0x4E && data[3]==0x47) {
		std::size_t pos = 8;
		while (pos + 12 <= sz) {
			const uint32_t chunkLen = (static_cast<uint32_t>(data[pos])   << 24) |
			                          (static_cast<uint32_t>(data[pos+1]) << 16) |
			                          (static_cast<uint32_t>(data[pos+2]) << 8)  |
			                           static_cast<uint32_t>(data[pos+3]);
			if (pos + 12 + chunkLen > sz) break;
			if (data[pos+4]=='p' && data[pos+5]=='H' && data[pos+6]=='Y' && data[pos+7]=='s'
			    && chunkLen == 9) {
				const uint32_t xdpm = (static_cast<uint32_t>(data[pos+8])  << 24) |
				                      (static_cast<uint32_t>(data[pos+9])  << 16) |
				                      (static_cast<uint32_t>(data[pos+10]) << 8)  |
				                       static_cast<uint32_t>(data[pos+11]);
				const uint8_t unit = data[pos + 16];
				if (unit == 1) return xdpm;
			}
			if (data[pos+4]=='I' && data[pos+5]=='D' && data[pos+6]=='A' && data[pos+7]=='T') break;
			pos += 12 + chunkLen;
		}
		return 0;
	}

	// ---- BMP (42 4D) ----
	if (data[0] == 0x42 && data[1] == 0x4D) {
		if (sz < 42) return 0;
		// The biXPelsPerMeter field only exists in BITMAPINFOHEADER (40 B)
		// and later headers (V2=52, V3=56, V4=108, V5=124). The legacy
		// BITMAPCOREHEADER is 12 B and has no DPI fields at all — bytes
		// 38..41 there belong to pixel data and must not be interpreted
		// as DPI.
		const uint32_t dibSize = static_cast<uint32_t>(data[14]) |
		                         (static_cast<uint32_t>(data[15]) << 8)  |
		                         (static_cast<uint32_t>(data[16]) << 16) |
		                         (static_cast<uint32_t>(data[17]) << 24);
		if (dibSize < 40) return 0;
		const int32_t xppm = static_cast<int32_t>(data[38]) |
		                     (static_cast<int32_t>(data[39]) << 8)  |
		                     (static_cast<int32_t>(data[40]) << 16) |
		                     (static_cast<int32_t>(data[41]) << 24);
		return xppm > 0 ? static_cast<unsigned int>(xppm) : 0;
	}

	// ---- TIFF (49 49 2A 00 / 4D 4D 00 2A) ----
	if ((data[0]==0x49 && data[1]==0x49 && data[2]==0x2A && data[3]==0x00) ||
	    (data[0]==0x4D && data[1]==0x4D && data[2]==0x00 && data[3]==0x2A)) {
		return ReadTiffDotsPerMeterX(data, 0, sz);
	}

	return 0;
}

} // namespace

void OPENLQM_API_IMPL OpenLQM::Fingerprint::LoadFromFilePath(const std::string& filePath, OpenLQM::PixelDensity resolutionOverride) {
	FILE* pFile = fopen(filePath.c_str(), "rb");
	if (!pFile) {
		throw std::invalid_argument(std::string("Failed to open {") + filePath + "}");
	}
	std::vector<unsigned char> fileBytes;
	char buf[2049];
	buf[2048] = 0;
	std::size_t numRead = 0;
	while ((numRead = fread(buf, 1, 2048, pFile)) > 0) {
		std::size_t lastSize = fileBytes.size();
		fileBytes.resize(numRead + lastSize);
		for (std::size_t i = 0; i < numRead; ++i) {
			fileBytes[i + lastSize] = static_cast<unsigned char>(buf[i]);
		}
	}
	const unsigned int dpmX = ReadDotsPerMeterX(fileBytes);
	const float INCHES_PER_METER = 39.3701f;
	fclose(pFile);
	float ppiF = OpenLQM::Core::ClampResolution(dpmX / INCHES_PER_METER);
	int ppi = std::lrint(ppiF);
	if (resolutionOverride != OpenLQM::PixelDensity::ppiInvalid) {
		ppi = static_cast<int>(static_cast<unsigned int>(resolutionOverride));
	}

	cv::Mat inputMat = cv::imdecode(fileBytes, cv::IMREAD_GRAYSCALE);
	if (inputMat.empty()) {
		throw std::invalid_argument(std::string("LoadFromFilePath failed to decode file (") + filePath + ")");
	}

	this->width = static_cast<unsigned int>(inputMat.cols);
	this->height = static_cast<unsigned int>(inputMat.rows);
	if (ppi == 500) {
		this->resolution = OpenLQM::PixelDensity::ppi500;
	}
	else if (ppi == 1000) {
		this->resolution = OpenLQM::PixelDensity::ppi1000;
	} else if (ppi == 2000) {
		this->resolution = OpenLQM::PixelDensity::ppi2000;
	} else if (ppi == 4000) {
		this->resolution = OpenLQM::PixelDensity::ppi4000;
	} else if (resolutionOverride != OpenLQM::PixelDensity::ppiInvalid) {
		this->resolution = resolutionOverride;
	} else {
		throw std::invalid_argument(std::string("Invalid input resolution (") + std::to_string(ppi) + "). Must be approximately 500, 1000, 2000, or 4000 pixels per inch");
	}
	this->buffer.resize(this->width * this->height);
	memcpy(this->buffer.data(), inputMat.data, this->buffer.size());

	this->bitDepth = PixelBitDepth::depth8;
	this->bitsPerPixel = 8;
}

void OpenLQM::Fingerprint::SetRoi(const std::vector<OpenLQM::Coordinate>& coordinates) {
	this->roi = coordinates;
}

void OpenLQM::BinarizeMap(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, unsigned char threshold) {
	dest.reserve(src.size());
	dest.resize(0);
	for (const unsigned char c : src) {
		dest.emplace_back(c >= threshold);
	}
}

void OpenLQM::AreaCellCount(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius) {
	dest.resize(src.size());

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			int index = y * width + x;

			if (src[static_cast<std::size_t>(index)] == 0) {
				dest[static_cast<std::size_t>(index)] = 0;
				continue;
			}

			int cellTotal = 0;
			int cellSum = 0;
			for (int y1 = y - radius; y1 <= y + radius; ++y1) {
				for (int x1 = x - radius; x1 <= x + radius; ++x1) {
					if (x1 >= 0 && x1 < width && y1 >= 0 && y1 < height) {
						cellSum += src[static_cast<std::size_t>(y1 * width + x1)];
						++cellTotal;
					}
				}
			}
			double val = 255.0 * static_cast<double>(cellSum) / static_cast<double>(cellTotal);
			// This rounding function is used to match the behavior of Visual Basic's
			//  conversion from floating point types to integral types (banker's rounding)
			dest[static_cast<std::size_t>(index)] = static_cast<unsigned char>(std::lrint(val));
		}
	}
}

void OpenLQM::DilateErodeCore(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius, unsigned char val) {
	dest = src;

	for (int y = 0; y < height; ++y) {
		int yClamped = std::clamp(y, 0, height - 1);
		int x = 0;
		while (x < width) {
			int xClamped = std::clamp(x, 0, width - 1);
			int clampedIndex = yClamped * width + xClamped;
			if (src[static_cast<std::size_t>(clampedIndex)] != val) {
				++x;
				continue;
			}

			for (int y1 = y - radius; y1 <= y + radius; ++y1) {
				int y1Clamped = std::clamp(y1, 0, height - 1);
				for (int x1 = x - radius; x1 <= x + radius; ++x1) {
					int x1Clamped = std::clamp(x1, 0, width - 1);
					int clampedInnerIndex = y1Clamped * width + x1Clamped;
					dest[static_cast<std::size_t>(clampedInnerIndex)] = val;
				}
			}
			x += 2; // Already marked the next row
		}
	}
}

void OpenLQM::Dilate(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius) {
	DilateErodeCore(src, dest, width, height, radius, 1);
}

void OpenLQM::Erode(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius) {
	DilateErodeCore(src, dest, width, height, radius, 0);
}

void OpenLQM::Open(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius) {
	//todo: refactor/replace these static scratch vectors
	static std::vector<unsigned char> scratchVec;

	Erode(src, scratchVec, width, height, radius);
	Dilate(scratchVec, dest, width, height, radius);
}

void OpenLQM::Close(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius) {
	static std::vector<unsigned char> scratchVec;

	Dilate(src, scratchVec, width, height, radius);
	Erode(scratchVec, dest, width, height, 1);
}

void OpenLQM::OpenClose(const std::vector<unsigned char>& src, std::vector<unsigned char>& dest, int width, int height, int radius) {
	static std::vector<unsigned char> scratchVec;

	Open(src, scratchVec, width, height, radius);
	Close(scratchVec, dest, width, height, 1);
}

int OpenLQM::CountNonZero(const std::vector<unsigned char>& src) {
	int count = 0;
	for (unsigned char val : src) {
		count += (val > 0);
	}
	return count;
}

void OpenLQM::Core::FloodFillImpl(std::vector<unsigned char>& src, int width, int height, int x, int y, unsigned char targetValue, unsigned char newValue) {
	if (x >= 0 && x < width && y >= 0 && y < height) {
		int index = y * width + x;
		if (src[static_cast<std::size_t>(index)] != targetValue) {
			return;
		}

		src[static_cast<std::size_t>(index)] = newValue;

		// Fill current row to left
		int left = x;
		int testIndex = y * width + left;
		while (left > 0 && src[static_cast<std::size_t>(testIndex - 1)] == targetValue) {
			testIndex -= 1;
			left -= 1;
			src[static_cast<std::size_t>(testIndex)] = newValue;
		}

		// Fill current row to right
		int right = x;
		testIndex = y * width + right;
		while (right < width - 1 && src[static_cast<std::size_t>(testIndex + 1)] == targetValue) {
			++testIndex;
			++right;
			src[static_cast<std::size_t>(testIndex)] = newValue;
		}

		// Check row above, recurse only if needed
		if (y > 0) {
			testIndex = (y - 1) * width + left;
			for (int testX = left; testX <= right; ++testX) {
				if (src[static_cast<std::size_t>(testIndex)] == targetValue) {
					FloodFillImpl(src, width, height, testX, y - 1, targetValue, newValue);
				}
				++testIndex;
			}
		}

		// Check row below, recurse only if needed
		if (y < height - 1) {
			testIndex = (y + 1) * width + left;
			for (int testX = left; testX <= right; ++testX) {
				if (src[static_cast<std::size_t>(testIndex)] == targetValue) {
					FloodFillImpl(src, width, height, testX, y + 1, targetValue, newValue);
				}
				++testIndex;
			}
		}
	}
}

void OpenLQM::Core::FloodFill(std::vector<unsigned char>& src, int width, int height, int x, int y, unsigned char targetValue, unsigned char newValue) {
	if (targetValue == newValue) {
		throw std::invalid_argument(std::string("FloodFill with identical values (") + std::to_string(targetValue) + ")");
	}

	FloodFillImpl(src, width, height, x, y, targetValue, newValue);
}

void OpenLQM::MarkRegions(std::vector<unsigned char>& src, int& regionCount, int width, int height) {
	// Note: we save the regions in a byte array- but that limits us to 255 regions

	for (int y = 0; y < height; ++y) {
		int rowIndex = y * width;
		for (int x = 0; x < width; ++x) {
			if (src[static_cast<std::size_t>(rowIndex + x)] != 1) {
				continue;
			}

			if ((regionCount + 1) >= 255) {
				regionCount = 255; // With 255 regions we givue up.With pixel noise this is not unreasonable
				return;
			}
			++regionCount;

			// Linear Recursion
			OpenLQM::Core::FloodFill(src, width, height, x, y, 1, static_cast<unsigned char>(regionCount));
		}
	}
}

void OpenLQM::RemoveRegions(std::vector<unsigned char>& src, int width, int height, const std::vector<bool>& keepRegions) {
	for (int y = 0; y < height; ++y) {
		int rowIndex = y * width;
		for (int x = 0; x < width; ++x) {
			int index = rowIndex + x;
			if (!keepRegions[src[static_cast<std::size_t>(index)]]) {
				src[static_cast<std::size_t>(index)] = 0;
			}
		}
	}
}

void OpenLQM::Mask(std::vector<unsigned char>& target, const std::vector<unsigned char>& mask) {
	auto targetIt = target.begin();
	auto maskIt = mask.begin();
	for (; targetIt != target.end(); ++targetIt, ++maskIt) {
		if (!(*maskIt)) {
			*targetIt = 0;
		}
	}
}

int OpenLQM::SumAllCells(const std::vector<unsigned char>& src) {
	int count = 0;
	for (unsigned char val : src) {
		count += val;
	}
	return count;
}

int OpenLQM::ScaleValues(int inValue, int inMin, int inMax, int outMin, int outMax) {
	// Linearly rescale a value in [inMin, outMin] to [outMin, outMax]
	double inDouble = static_cast<double>(inValue - inMin) / static_cast<double>(inMax - inMin);
	double scaled = inDouble * static_cast<double>(outMax - outMin) + static_cast<double>(outMin);
	double clamped = std::clamp(scaled, static_cast<double>(outMin), static_cast<double>(outMax));

	// This rounding function is used to match the behavior of Visual Basic's
	//  conversion from floating point types to integral types (banker's rounding)
	return std::lrint(clamped);
}

OpenLQM::Core::ImageAttributes::ImageAttributes() {

}

OpenLQM::Core::ImageAttributes::ImageAttributes(const std::vector<unsigned char>& buf, int width, int height, int resolution) {
	Load(buf, width, height, resolution);
}

void OpenLQM::Core::ImageAttributes::Init(const std::vector<unsigned char>& buf, int width, int height, int resolution) {
	*this = {};
	Load(buf, width, height, resolution);
}

void OpenLQM::Core::ImageAttributes::Load(const std::vector<unsigned char>& buf, int width, int height, int resolution) {
	Width = width;
	Height = height;
	Resolution = resolution;
	PixelCount = width * height;
	Recalc(buf);
}

void OpenLQM::Core::ImageAttributes::CalcStatisticsFromHistogram() {
	int64_t RunningTotal = 0, RunningProduct = 0;
	int QTargets[5] = {};

	for (int Qcount = 0; Qcount <= 4; ++Qcount) {
		Quartile[Qcount] = -1;
		QTargets[Qcount] = static_cast<int>(std::lrint(static_cast<double>(Qcount * PixelCount) / 4.0));
	}
	QTargets[0] = 1; // to compare to RunningTotal, either min's target needs to be 1 with a comarison of >=, or max's target needs to be PixelCount-1 with a comparison of >

	for (int valBin = 0; valBin < static_cast<int>(Hist.size()); ++valBin) {
		int valCount = Hist[static_cast<std::size_t>(valBin)];
		if (valCount > 0) {
			++GrayValueCount;
		}

		RunningTotal += valCount;
		RunningProduct += static_cast<int64_t>(valCount) * valBin;
		for (int Qcount = 0; Qcount <= 4; ++Qcount) {
			if (Quartile[Qcount] < 0 && RunningTotal >= QTargets[Qcount]) {
				Quartile[Qcount] = valBin;
			}
		}
	}

	Range = Quartile[4] - Quartile[0];
	Mean = static_cast<int>(std::lrint(static_cast<double>(RunningProduct) / static_cast<double>(PixelCount)));
}

void OpenLQM::Core::ImageAttributes::Recalc(const std::vector<unsigned char>& buf) {
	for (int& i : Hist) {
		i = 0;
	}
	PixelCount = static_cast<int>(buf.size());

	for (unsigned char val : buf) {
		Hist[val] += 1;
	}
	CalcStatisticsFromHistogram();
}

OpenLQM::Core::AggregateQuality::AggregateQuality(const std::vector<unsigned char>& LocQMap, int width, int height, OpenLQM::Core::ImageAttributes& ImgAttrib) {
	OpenLQM::Core::ImageAttributes LocQAttribs(LocQMap, width, height, 0);
	for (int Lcounter = 0; Lcounter < RESULT_VECTOR_LENGTH; ++Lcounter) {
		for (int counter = Lcounter; counter < RESULT_VECTOR_LENGTH; ++counter) {
			m_totalArea[static_cast<std::size_t>(Lcounter)] += LocQAttribs.Hist[static_cast<std::size_t>(counter)]; // cumulative
		}
	}

	// For the below, values for 0 and 1 are ignored
	std::array<std::vector<unsigned char>, RESULT_VECTOR_LENGTH> LocQBin; // binarized
	std::array<std::vector<unsigned char>, RESULT_VECTOR_LENGTH> LocQEntropy; // "Entropy" - proportion of neighboring nonzero cells
	std::array<std::vector<unsigned char>, RESULT_VECTOR_LENGTH> LocQOC; // OpenClose - gets rid of minor holes, small areas, indents, protrusions
	std::array<std::vector<unsigned char>, RESULT_VECTOR_LENGTH> LocQregions;
	std::array<int, RESULT_VECTOR_LENGTH> LCA; // the index # of the largest area in LocQregions
	std::array<OpenLQM::Core::ImageAttributes, RESULT_VECTOR_LENGTH> LocQregionAttribs;
	std::vector<bool> keepRegions; // temp flag for regions to be retained (LCA or over a given size)
	std::array<std::vector<unsigned char>, RESULT_VECTOR_LENGTH> ConsistencyGFA;
	std::array<std::vector<unsigned char>, RESULT_VECTOR_LENGTH> ConsistencyLCA;

	// No further stats for LocQ 0 & 1
	//int regionCount = 1;
	for (int Lcounter = 0; Lcounter < RESULT_VECTOR_LENGTH; ++Lcounter) {
		// Binarize (by local quality)
		OpenLQM::BinarizeMap(LocQMap, LocQBin[static_cast<std::size_t>(Lcounter)], static_cast<unsigned char>(Lcounter));

		// "entropy"
		// For a box with the given "radius" (e.g. rad=8 means a 17x17 mapbox=68x68@500ppi=~3ridges in each direction),
		//     indicate the proportion that are nonzero on a 0-255 scale: 0=no nonzero cells..255=all nonzero cells
		OpenLQM::AreaCellCount(LocQBin[static_cast<std::size_t>(Lcounter)], LocQEntropy[static_cast<std::size_t>(Lcounter)], width, height, 8);

		// Filter (open close) with a "radius"=1 or 2 -3x3 or 5x5 box
		OpenLQM::OpenClose(LocQBin[static_cast<std::size_t>(Lcounter)], LocQOC[static_cast<std::size_t>(Lcounter)], width, height, m_openCloseRadii[static_cast<std::size_t>(Lcounter)]);
		m_totalFilteredArea[static_cast<std::size_t>(Lcounter)] = OpenLQM::CountNonZero(LocQOC[static_cast<std::size_t>(Lcounter)]);

		// Mark regions: all contiguous areas are marked with the same non-zero number
		LocQregions[static_cast<std::size_t>(Lcounter)] = LocQOC[static_cast<std::size_t>(Lcounter)];
		int regionCount = 1;
		OpenLQM::MarkRegions(LocQregions[static_cast<std::size_t>(Lcounter)], regionCount, width, height); // RECURSIVE - can fail with stack size problems
		m_regionCount[static_cast<std::size_t>(Lcounter)] = regionCount;

		// Mask out small regions
		LCA[static_cast<std::size_t>(Lcounter)] = 0;
		keepRegions.resize(static_cast<std::size_t>(m_regionCount[static_cast<std::size_t>(Lcounter)] + 2)); // Ignore values 0, 1
		LocQregionAttribs[static_cast<std::size_t>(Lcounter)] = OpenLQM::Core::ImageAttributes(LocQregions[static_cast<std::size_t>(Lcounter)], width, height);
		for (int counter = 1; counter <= m_regionCount[static_cast<std::size_t>(Lcounter)]; ++counter) {
			keepRegions[static_cast<std::size_t>(counter + 1)] = (LocQregionAttribs[static_cast<std::size_t>(Lcounter)].Hist[static_cast<std::size_t>(counter + 1)] > m_ignoreThreshold[static_cast<std::size_t>(Lcounter)]);
			if (LocQregionAttribs[static_cast<std::size_t>(Lcounter)].Hist[static_cast<std::size_t>(counter + 1)] > m_largestContiguousArea[static_cast<std::size_t>(Lcounter)]) {
				LCA[static_cast<std::size_t>(Lcounter)] = counter + 1;
				m_largestContiguousArea[static_cast<std::size_t>(Lcounter)] = LocQregionAttribs[static_cast<std::size_t>(Lcounter)].Hist[static_cast<std::size_t>(counter) + 1];
			}
		}

		if (Lcounter == 2) {
			keepRegions[static_cast<std::size_t>(LCA[static_cast<std::size_t>(Lcounter)])] = true; // Keep the largest region JUST for L2 AND all regions bigger than m_ignoreThreshold (all levels)
		}
		OpenLQM::RemoveRegions(LocQregions[static_cast<std::size_t>(Lcounter)], width, height, keepRegions);
		LocQregionAttribs[static_cast<std::size_t>(Lcounter)].Recalc(LocQregions[static_cast<std::size_t>(Lcounter)]);
	}

	// Now what we care about is level 3&4 in a) GFA and b) LCA

    // a) Mask GFA
    //       Mask the original image to ignore everything that isn't an acceptable L2 region
    // GFA=Good flow areas, means large areas of L2 or better
	GoodFlowAreaOfRidgeFlow = std::vector<unsigned char>(LocQMap.begin(), LocQMap.end());
	OpenLQM::Mask(GoodFlowAreaOfRidgeFlow, LocQregions[2]); // Mask everything that isn't an acceptable L2 region
    OpenLQM::Core::ImageAttributes LocQAttribs_GFA(GoodFlowAreaOfRidgeFlow, width, height);

	// b) Mask LCA (a subset of GFA for complicated images)
	std::vector<unsigned char> LocQMap_LCAMasked = std::vector<unsigned char>(LocQMap.begin(), LocQMap.end());
	keepRegions.resize(static_cast<std::size_t>(m_regionCount[2] + 2)); // Ignore values 0, 1
	for (int counter = 1; counter <= (static_cast<int>(keepRegions.size())) - 1; ++counter) {
		keepRegions[static_cast<std::size_t>(counter)] = (counter == LCA[2]);
	}

	LargestContiguousAreaOfRidgeFlow = std::vector<unsigned char>(LocQregions[2].begin(), LocQregions[2].end());
	OpenLQM::RemoveRegions(LargestContiguousAreaOfRidgeFlow, width, height, keepRegions);
	OpenLQM::Mask(LocQMap_LCAMasked, LargestContiguousAreaOfRidgeFlow); // Mask everything that isn't an acceptable L2 region
	OpenLQM::Core::ImageAttributes LocQAttribs_LCA(LocQMap_LCAMasked, width, height);

	for (int Lcounter = 0; Lcounter < RESULT_VECTOR_LENGTH; ++Lcounter) {
		for (int counter = Lcounter; counter < RESULT_VECTOR_LENGTH; ++counter) {
			m_areaWithinGFA[static_cast<std::size_t>(Lcounter)] += LocQAttribs_GFA.Hist[static_cast<std::size_t>(counter)]; // cumulative
			m_areaWithinLCA[static_cast<std::size_t>(Lcounter)] += LocQAttribs_LCA.Hist[static_cast<std::size_t>(counter)]; // cumulative
		}

		// Take the entropy maps (indicating consistency in a 0..255 map), mask using L2 regions, and sum
		// - results in a count of how much L3 and L4 is in a consistent region
		ConsistencyGFA[static_cast<std::size_t>(Lcounter)] = LocQEntropy[static_cast<std::size_t>(Lcounter)];
		OpenLQM::Mask(ConsistencyGFA[static_cast<std::size_t>(Lcounter)], LocQregions[2]); // Mask everything that isn't GFA

		ConsistencyLCA[static_cast<std::size_t>(Lcounter)] = LocQEntropy[static_cast<std::size_t>(Lcounter)];
		OpenLQM::Mask(ConsistencyLCA[static_cast<std::size_t>(Lcounter)], LargestContiguousAreaOfRidgeFlow); // Mask everything that isn't LCA

		// ConsistencyGFA is a buffer with each cell 0..255
		// m_consistencyWithinGFA is an integer counting pixels, devaluing pixels with neighbors outside region
		int sumI = OpenLQM::SumAllCells(ConsistencyGFA[static_cast<std::size_t>(Lcounter)]);
		long long sumLL = static_cast<long long>(sumI);
		double sumSS = static_cast<double>(sumLL);
		sumSS = sumSS / 255.0;
		m_consistencyWithinGFA[static_cast<std::size_t>(Lcounter)] = static_cast<int>(std::lrintl(sumSS)); // Dividing by 255 makes the units pixels again
		sumI = OpenLQM::SumAllCells(ConsistencyLCA[static_cast<std::size_t>(Lcounter)]);
		sumLL = sumI;
		sumSS = static_cast<double>(sumLL);
		//m_consistencyWithinLCA[Lcounter] = (int)std::lrintl(sum / 255.0); // Dividing by 255 makes the units pixels again
		m_consistencyWithinLCA[static_cast<std::size_t>(Lcounter)] = static_cast<int>(sumSS / 255.0); // Dividing by 255 makes the units pixels again

	}

	CalcOverallQuality_Consist();

	if (ImgAttrib.Range == 0) {
		OverallClarity = 0;
	}
	else {
		OverallClarity = (OverallQuality_GFA_Consist + OverallQuality_LCA_Consist) / 2;
	}
}

void OpenLQM::Core::AggregateQuality::CalcOverallQuality_Consist() {
	OverallQuality_GFA_Consist = CalcOverallQuality_Consist(m_consistencyWithinGFA);
	OverallQuality_LCA_Consist = CalcOverallQuality_Consist(m_consistencyWithinLCA);
}

int OpenLQM::Core::AggregateQuality::CalcOverallQuality_Consist(std::array<int, OpenLQM::Core::RESULT_VECTOR_LENGTH>& ConsistencyArray) {
	// Note: This algorithm combines the L3 (minutiae confidence) and L4 (ridge edge) in a weighted average
	double ConsistencyArray34 = (2.0 * ConsistencyArray[3] + ConsistencyArray[4]) / 3.0; // weighted average of 2xL3 + L4, rescaled to units of pix@125ppi
	// Note: These odd-seeming numbers are fairly normal-looking in sq.in. (si)
	if (m_totalArea[1] == 0) {
		return 1; // 1: apparently blank; 0: reserved for completely blank
	}
	else if (ConsistencyArray[2] == 0) {
		return OpenLQM::ScaleValues(m_totalArea[1], 0, 10937, 2, 9); // :unusable: 10937 = 0.7si
	}
	else if (ConsistencyArray34 <= 156.0) { // 0.01si
		double inValue = static_cast<double>(ConsistencyArray[2]) + ConsistencyArray34;
		return OpenLQM::ScaleValues(std::lrint(inValue), 0, 10937, 10, 19); // 10-19:exclusion only: 10937 = 0.7si - but doubleweight 34
	}
	else if (ConsistencyArray34 <= 1562.0) { // 0.1si
		return OpenLQM::ScaleValues(std::lrint(ConsistencyArray34), 156, 1562, 20, 40); // Difficult - note overlap at 40
	}
	else if (ConsistencyArray34 <= 3125.0) { // 0.2si
		return OpenLQM::ScaleValues(std::lrint(ConsistencyArray34), 1562, 3125, 40, 60); // Moderate - note overlap at 60
	}
	else if (ConsistencyArray34 <= 6250.0) { // 0.4si
		return OpenLQM::ScaleValues(std::lrint(ConsistencyArray34), 3125, 6250, 60, 80); // Straightforward - note overlap at 80
	}
	else if (ConsistencyArray34 <= 10937.0) { // 0.7si
		return OpenLQM::ScaleValues(std::lrint(ConsistencyArray34), 6250, 10937, 80, 90); // Easy - note overlap at 90
	}
	else { // Ideal
		if (ConsistencyArray[4] <= 6250) { // 0.4si - no values >=95 without L4
			return OpenLQM::ScaleValues(std::lrint(ConsistencyArray34), 10937, 23437, 90, 95);
		}
		else {
			return OpenLQM::ScaleValues(std::lrint(ConsistencyArray34), 10937, 23437, 90, 99); // Notice the overlap in 90-95 range
		}
	}
}

bool OpenLQM::Core::CheckHeuristics(double LCAQ1plus, double LCAQ3plus, int minQ2plus) {
	return !(LCAQ1plus < 30.0 || LCAQ3plus <= 0.5 || minQ2plus > 795);
}

double OpenLQM::Core::CalculateLQMetric12(double predictedMatchScore) {
	if (predictedMatchScore >= 20183.3635180282) {
		return 1.0;
	}
	else if (predictedMatchScore >= 17777.6337566414) {
		return predictedMatchScore * 0.00001731976 + 0.650428949;
	}
	else if (predictedMatchScore >= 16349.8836341313) {
		return predictedMatchScore * 0.00000767985 + 0.821803704;
	}
	else if (predictedMatchScore >= 14199.1879480481) {
		return predictedMatchScore * 0.000039908 + 0.294877295;
	}
	else if (predictedMatchScore >= 14038.5648550078) {
		return predictedMatchScore * 0.00038312337 + -4.57850235;
	}
	else if (predictedMatchScore >= 12822.8383519244) {
		return predictedMatchScore * 0.00002827527 + 0.403055732;
	}
	else if (predictedMatchScore >= 11461.0338401133) {
		return (predictedMatchScore * 0.0000474327) + 0.157403087;
	}
	else if (predictedMatchScore >= 10592.4798643788) {
		return (predictedMatchScore * 0.00006214019) + -0.01115983;
	}
	else if (predictedMatchScore >= 9830.0511721355) {
		return (predictedMatchScore * 0.00013187751) + -0.74985102;
	}
	else if (predictedMatchScore >= 9178.87500654071) {
		return (predictedMatchScore * 0.0000714271) + -0.155620463;
	}
	else if (predictedMatchScore >= 9070.83210405052) {
		return (predictedMatchScore * 0.00027222301) + -1.998701032;
	}
	else if (predictedMatchScore >= 8790.52152277175) {
		return (predictedMatchScore * 0.00010912266) + -0.519245061;
	}
	else if (predictedMatchScore >= 8075.4871812788) {
		return (predictedMatchScore * 0.00005812658) + -0.070962942;
	}
	else if (predictedMatchScore >= 8060.08296383831) {
		return (predictedMatchScore * 0.00422638585) + -33.73168728;
	}
	else if (predictedMatchScore >= 7858.24815900676) {
		return (predictedMatchScore * 0.00037748929) + -2.709261623;
	}
	else if (predictedMatchScore >= 7497.88248290373) {
		return (predictedMatchScore * 0.0001422458) + -0.860659946;
	}
	else if (predictedMatchScore >= 7349.42804513656) {
		return (predictedMatchScore * 0.00042454245) + -2.977287014;
	}
	else if (predictedMatchScore >= 7120.75746589422) {
		return (predictedMatchScore * 0.00002839677) + -0.065842896;
	}
	else if (predictedMatchScore >= 6789.54131448483) {
		return (predictedMatchScore * 0.00008822269) + -0.491848729;
	}
	else if (predictedMatchScore >= 6546.9798744552) {
		return (predictedMatchScore * 0.00002944762) + -0.092792976;
	}
	else if (predictedMatchScore >= 5409.8752849724) {
		return (predictedMatchScore * 0.00001147078) + 0.024901025;
	}
	else if (predictedMatchScore >= 2849.35331665492) {
		return (predictedMatchScore * 0.00003396047) + -0.096765369;
	}
	else if (predictedMatchScore >= 2849.35331665492 && predictedMatchScore < 5409.8752849724) {
		return (predictedMatchScore * 0.00003396047) + -0.096765369;
	}
	else {
		return 0.0;
	}
}

int OpenLQM::Core::PredictValueDeterminations(double predictedMatchScore, bool determFlag) {
	static std::vector<std::tuple<double, double, double>> MATCH_VECTORS_TRUE = {
		std::tuple<double, double, double>(14826.3514217431, 0.0, 1.0),
		std::tuple<double, double, double>(13311.3968181226, 0.00000904227, 0.865936112),
		std::tuple<double, double, double>(11316.9853097232, 0.00000411108, 0.931577187),
		std::tuple<double, double, double>(10748.4893430339, 0.00002320141, 0.715532152),
		std::tuple<double, double, double>(10002.5948865749, 0.00000084002, 0.955883332),
		std::tuple<double, double, double>(9813.5219210502, 0.00010493975, -0.085384098),
		std::tuple<double, double, double>(9758.4375301615, 0.00050427675, -4.004286525),
		std::tuple<double, double, double>(9151.30021647409, 0.00001960802, 0.725323007),
		std::tuple<double, double, double>(8896.03487570362, 0.00000909987, 0.8214863),
		std::tuple<double, double, double>(8488.62592428182, 0.00012539654, -0.213092972),
		std::tuple<double, double, double>(8301.97874339333, 0.00009653517, 0.031900374),
		std::tuple<double, double, double>(8259.77851100455, 0.00078988507, -5.724275692),
		std::tuple<double, double, double>(7982.38590374231, 0.00003353487, 0.523009385),
		std::tuple<double, double, double>(7809.79906551272, 0.0001919145, -0.741237929),
		std::tuple<double, double, double>(7672.98918017373, 0.00074850837, -5.088124183),
		std::tuple<double, double, double>(7421.11243444695, 0.00005616546, 0.224215463),
		std::tuple<double, double, double>(7185.86012982543, 0.00020095719, -0.850300227),
		std::tuple<double, double, double>(6890.45439098391, 0.00002553472, 0.410261059),
		std::tuple<double, double, double>(4971.91832907826, 0.00000770292, 0.533130294),
		std::tuple<double, double, double>(4084.19262471953, 0.00008046244, 0.171375918),
		std::tuple<double, double, double>(2849.35331665492, 0.00040491098, -1.153734457)
	};

	static std::vector<std::tuple<double, double, double>> MATCH_VECTORS_FALSE = {
		std::tuple<double, double, double>(12982.5296597007, 0.0, 1.0),
		std::tuple<double, double, double>(11316.9853097232, 0.00000504541, 0.934497759),
		std::tuple<double, double, double>(8896.03487570362, 0.00000304406, 0.957147099),
		std::tuple<double, double, double>(8442.61676793259, 0.00001650341, 0.837412251),
		std::tuple<double, double, double>(7806.50006743808, 0.00003300019, 0.698136206),
		std::tuple<double, double, double>(7666.15857084569, 0.00019367603, -0.556179724),
		std::tuple<double, double, double>(4971.91832907826, 0.00000507669, 0.889652749),
		std::tuple<double, double, double>(3584.46167254856, 0.00024755011, -0.315905304),
		std::tuple<double, double, double>(2849.35331665492, 0.00077733924, -2.214914145)
	};

	auto& predictionArray = (determFlag ? MATCH_VECTORS_TRUE : MATCH_VECTORS_FALSE);
	std::tuple<double, double, double> selectedVec{0.0, 0.0, 0.0};

	for (const std::tuple<double, double, double>& predVec : predictionArray) {
		if (predictedMatchScore >= std::get<0>(predVec)) {
			selectedVec = predVec;
			break;
		}
	}

	double predictedScore = (predictedMatchScore * std::get<1>(selectedVec) + std::get<2>(selectedVec)) * 100.0;
	return static_cast<int>(std::lrint(predictedScore));

}

long OpenLQM::Core::EmbedMinutiaeDetail(
	const LQMParams* pLqmParams,
	unsigned char* grayImage,
	unsigned char* binaryImage,
	int inputResolution,
	int imgWidth, int imgHeight,
	char* DirMap,
	char* QualMap,
	int maxMinutiae,
	MinutiaOut* MinutiaeArray,
	unsigned char* gray10pctMap,  // grayscale of 10th %ile (dark)
	unsigned char* gray90pctMap,  // grayscale of 90th %ile (light)
	unsigned char* grayMedianMap, // median grayscale
	unsigned char* grayRangeMap, // median grayscale
	unsigned char* grayCountMap,   // number of gray values in use
	double* maxMagnitude,
	double* normMagnitude,
	double* lowFreqMagnitude,
	unsigned char* directionChangeMap, // changes in direction map
	unsigned char* validNeighbors,
	unsigned char* curvatureMap
) {
	int bw, bh;
	std::vector<int> direction_map, low_contrast_map, low_flow_map, high_curve_map;
	std::vector<unsigned char> binData;
	int map_w, map_h;

	if (pLqmParams == nullptr) {
		pLqmParams = &OpenLQM::Presets::LQM_PARAMS;
	}

	Minutiae minutiae;
	DetectMinutiae(
		minutiae,
		direction_map, low_contrast_map,
		low_flow_map, high_curve_map,
		&map_w, &map_h,
		binData, &bw, &bh,
		grayImage, imgWidth, imgHeight, *pLqmParams,
		gray10pctMap, gray90pctMap, grayMedianMap, grayRangeMap, grayCountMap,
		maxMagnitude, normMagnitude, lowFreqMagnitude,
		directionChangeMap, validNeighbors, curvatureMap
	);

	memcpy(binaryImage, binData.data(), static_cast<std::size_t>(bw*bh));

	for (int thisY = 0; thisY < map_h; ++thisY) {
		for (int thisX = 0; thisX < map_w; ++thisX) {
			int arrayPos = thisY*map_w+thisX;
			DirMap[static_cast<std::size_t>(arrayPos)] = static_cast<char>(direction_map[static_cast<std::size_t>(arrayPos)]);
			if (!QualMap) {
				continue;
			}

			if (low_contrast_map[static_cast<std::size_t>(arrayPos)] || direction_map[static_cast<std::size_t>(arrayPos)]<0) {
				QualMap[static_cast<std::size_t>(arrayPos)]=0;
			}
			else {
				//set baseline quality before looking at neighbors (will subtract QualOffset below)
				if (low_flow_map[static_cast<std::size_t>(arrayPos)] || high_curve_map[static_cast<std::size_t>(arrayPos)]) {
					QualMap[static_cast<std::size_t>(arrayPos)]=3; //offset will be -1..-2
				}
				else {
					QualMap[static_cast<std::size_t>(arrayPos)]=4; //offset will be 0..-2
				}

				if (thisY < LQM_NEIGHBOR_DELTA || thisY > map_h - 1 - LQM_NEIGHBOR_DELTA ||
					thisX < LQM_NEIGHBOR_DELTA || thisX > map_w - 1 - LQM_NEIGHBOR_DELTA)
				{
					QualMap[static_cast<std::size_t>(arrayPos)]=1;
				}
				else {
					int QualOffset = 0;
					for (int compY = thisY - LQM_NEIGHBOR_DELTA; compY <= thisY + LQM_NEIGHBOR_DELTA; ++compY) {
						for (int compX = thisX - LQM_NEIGHBOR_DELTA; compX <= thisX + LQM_NEIGHBOR_DELTA; ++compX) {
							int arrayPos2 = compY*map_w+compX;
							if (low_contrast_map[static_cast<std::size_t>(arrayPos2)] || direction_map[static_cast<std::size_t>(arrayPos2)] < 0) {
								QualOffset = -2;
								break;
							} else if (low_flow_map[static_cast<std::size_t>(arrayPos2)] || high_curve_map[static_cast<std::size_t>(arrayPos2)]) {
								QualOffset = std::min(QualOffset,-1);
							}
						}
					}
					QualMap[static_cast<std::size_t>(arrayPos)] += static_cast<char>(QualOffset);
				}
			}
		}
	}

	int ret = FilterAndOutputMinutiae(MinutiaeArray, maxMinutiae, minutiae, inputResolution);

	return ret;
}

int OpenLQM::Core::ProcessClarityMap(
	unsigned char* pOutClarityMapNonLCA, int mapWidth, int mapHeight,
	const unsigned char* pLocQ,
	const unsigned char* pGrayRangeMap,
	const unsigned char* pCurvatureMap,
	const unsigned char* pDirChangeMap,
	const double* pLowFreqMagMap,
	const double* pMaxMagMap,
	const unsigned char* pGrayMedianMap,
	const unsigned char* pGrayCountMap,
	const unsigned char* pValidNeighbors,
	const double* pNormMagMap,
	const unsigned char* pQualMap,
	const std::vector<OpenLQM::Coordinate>* roi
) {
	cv::Mat clarityMap(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(static_cast<const void*>(pLocQ)), static_cast<std::size_t>(mapWidth));
	cv::Mat m_clarityMapNonLCA(mapHeight, mapWidth, CV_8UC1, const_cast<void*>(static_cast<const void*>(pOutClarityMapNonLCA)), static_cast<std::size_t>(mapWidth));

	GenerateNewClarityMapRandomForest(
		mapWidth, mapHeight,
		clarityMap,
		pGrayRangeMap,
		pCurvatureMap,
		pDirChangeMap,
		pLowFreqMagMap,
		pMaxMagMap,
		pGrayMedianMap,
		pGrayCountMap,
		pValidNeighbors,
		pNormMagMap,
		pQualMap
	);

	clarityMap.copyTo(m_clarityMapNonLCA);
    cv::blur(m_clarityMapNonLCA, m_clarityMapNonLCA, cv::Size(3,3));

	if (roi && !roi->empty()) {
		CropToRoi(clarityMap, *roi);
	} else {
		CropToLargestContiguousRegionOfRidgeFlow(clarityMap);
	}

	return 0;
}

bool OpenLQM::Core::CalculateMetrics(
	const unsigned char* pInputImg, int inputWidth, int inputHeight, int resolution,
	const unsigned char* pLocQ, int mapWidth, int mapHeight,
	MinutiaOut* pMinutiae, int minutiaCount,
	OpenLQM::Metrics& metrics,
	OpenLQM::Supplement::QualityMeasures* pQualMeasures
) {
	std::vector<unsigned char> inputVec(pInputImg, pInputImg + (inputWidth * inputHeight));
	OpenLQM::Core::ImageAttributes imgAttributes(inputVec, inputWidth, inputHeight, resolution);

	std::vector<unsigned char> locqVec(pLocQ, pLocQ + (mapWidth * mapHeight));

	OpenLQM::Core::AggregateQuality aggQual(locqVec, mapWidth, mapHeight, imgAttributes);

	double AreaQ3plus = static_cast<double>(aggQual.m_totalArea[3]) / OpenLQM::Core::PIXEL_MAPPING;
	double AreaQ3plusDiv2 = static_cast<double>(aggQual.m_totalArea[3]) / (static_cast<double>(aggQual.m_totalArea[2] - aggQual.m_totalArea[3])) / OpenLQM::Core::PIXEL_MAPPING;
	double GFAQ3plus = static_cast<double>(aggQual.m_areaWithinGFA[3]) / OpenLQM::Core::PIXEL_MAPPING;
	double LCAQ1plus = static_cast<double>(aggQual.m_largestContiguousArea[1]) / OpenLQM::Core::AggregateQuality::MAP_TO_SQ_MM_FACTOR;
	double LCAQ3plus = static_cast<double>(aggQual.m_largestContiguousArea[3]) / OpenLQM::Core::AggregateQuality::MAP_TO_SQ_MM_FACTOR;

	std::vector<MinutiaOut> sortedMinutiae;
	FilterMinutiaeByQuality(pLocQ, mapWidth, mapHeight, inputWidth, inputHeight, pMinutiae, minutiaCount, sortedMinutiae);
	std::array<int, 3> minutiaCounts;
	CalculateMinutiaCounts(sortedMinutiae, minutiaCounts);
	int minQ2plus = minutiaCounts[0] + minutiaCounts[1] + minutiaCounts[2];

	bool hFlag = OpenLQM::Core::CheckHeuristics(LCAQ1plus, LCAQ3plus, minQ2plus);

	int minQ2 = minutiaCounts[0];
	int minQ4 = minutiaCounts[2];

	if (std::isnan(AreaQ3plusDiv2)) {
		AreaQ3plusDiv2 = 0.0;
	}

	if (std::isinf(AreaQ3plusDiv2)) {
		AreaQ3plusDiv2 = 10.89;
	}

	double predictedMatchScore = 7666.16 + (AreaQ3plusDiv2 * 32.38) + (AreaQ3plus * 17.71) + (GFAQ3plus * 20.87) - (static_cast<double>(minQ2) * 21.05) + (static_cast<double>(minQ4) * 7.71);

	int lqMetric12 = static_cast<int>(std::lrint(OpenLQM::Core::CalculateLQMetric12(predictedMatchScore) * 100.0));
	int lqm01Score = std::min(lqMetric12, 99);
	int vid = OpenLQM::Core::PredictValueDeterminations(predictedMatchScore, true);
	int vcmp = OpenLQM::Core::PredictValueDeterminations(predictedMatchScore, false);

	if (!hFlag) {
		lqm01Score = 0;
		vid = 0;
		vcmp = 0;
	}

	metrics.overallQuality = lqm01Score;
	metrics.vid = vid;
	metrics.vcmp = vcmp;
	metrics.overallClarity = aggQual.OverallClarity;
	metrics.areaOfImpression = aggQual.m_totalArea[1] / OpenLQM::Core::PIXEL_MAPPING;
	metrics.areaOfRidgeFlow = aggQual.m_totalArea[2] / OpenLQM::Core::PIXEL_MAPPING;
	metrics.areaOfGoodRidgeFlow = aggQual.m_totalArea[3] / OpenLQM::Core::PIXEL_MAPPING;
	metrics.areaOfClearLevel3Detail = aggQual.m_totalArea[4] / OpenLQM::Core::PIXEL_MAPPING;
	metrics.largestContiguousAreaOfRidgeFlow = aggQual.m_largestContiguousArea[2] / OpenLQM::Core::PIXEL_MAPPING;
	metrics.largestContiguousAreaOfGoodRidgeFlow = aggQual.m_largestContiguousArea[3] / OpenLQM::Core::PIXEL_MAPPING;
	metrics.automatedMinutiaeYellow = minutiaCounts[0];
	metrics.automatedMinutiaeGreen = minutiaCounts[1];
	metrics.automatedMinutiaeBlue = minutiaCounts[2];

	if (pQualMeasures) {
		pQualMeasures->Area3plusDiv2 = AreaQ3plusDiv2;
		pQualMeasures->Area3plus = AreaQ3plus;
		pQualMeasures->GFA3plus = GFAQ3plus;
		pQualMeasures->LCA1plus = LCAQ1plus;
		pQualMeasures->LCA3plus = LCAQ3plus;
		pQualMeasures->Min2 = minQ2;
		pQualMeasures->Min2plus = minQ2plus;
		pQualMeasures->Min4 = minQ4;
	}

	return hFlag;
}

void OPENLQM_API_IMPL OpenLQM::Supplement::PrintMetrics(const OpenLQM::Metrics& metrics, const std::string& inputPath, bool verbose, bool printHeaders) {
	std::vector<std::string> headers{
		"OQ",
		"VID",
		"VCMP",
		"OC",
		"Area_RYGB",
		"Area_YGB",
		"Area_GB",
		"Area_B",
		"LCA_YGB",
		"LCA_GRF",
		"AM_Y",
		"AM_G",
		"AM_B"
	};
	std::stringstream ss;
	ss << metrics.overallQuality << "  " << metrics.vid << "  " << metrics.vcmp << "  " << metrics.overallClarity << "  " << std::fixed << std::setprecision(1) << ROUND_OUT(metrics.areaOfImpression) << "  " << ROUND_OUT(metrics.areaOfRidgeFlow) << "  " << ROUND_OUT(metrics.areaOfGoodRidgeFlow) << "  " << ROUND_OUT(metrics.areaOfClearLevel3Detail) << "  " << ROUND_OUT(metrics.largestContiguousAreaOfRidgeFlow) << "  " << ROUND_OUT(metrics.largestContiguousAreaOfGoodRidgeFlow) << "  " << metrics.automatedMinutiaeYellow << "  " << metrics.automatedMinutiaeGreen << "  " << metrics.automatedMinutiaeBlue << std::endl;
	std::vector<std::string> scores;
	std::string score;
	std::string outStr = ss.str();
	std::stringstream tokenizer(outStr);
	while (tokenizer >> score) {
		if (!score.empty()) {
			scores.push_back(score);
		}
	}

	//verbose indices:
	int verboseCutoff = verbose ? static_cast<int>(headers.size()) : 3;

	if (printHeaders) {
		std::cout << "Filename\t";
		for (int i = 0; i < static_cast<int>(headers.size()) && i < verboseCutoff; ++i) {
			std::cout << headers[static_cast<std::size_t>(i)] << "\t";
		}
		std::cout << std::endl;
	}

	std::cout << inputPath << "\t";
	for (int i = 0; i < static_cast<int>(scores.size()) && i < verboseCutoff; ++i) {
		std::string sc = scores[static_cast<std::size_t>(i)];
		std::cout << sc << "\t";
	}
	std::cout << std::endl;
}

bool OpenLQM::GetAllMetricsAndQualityMapFromFingerprint_Base(const OpenLQM::Fingerprint& inpImg, OpenLQM::Metrics& metrics, OpenLQM::QualityMap& qualMap, bool inputQualityMap, bool justQualityMap, OpenLQM::Supplement::QualityMeasures* pQualMeasures, OpenLQM::Supplement::FeatureMaps* pFeatureMaps) {
	OpenLQM::Fingerprint normImg = OpenLQM::Core::CreateNormFingerprint(inpImg);

	cv::Mat grayImage = cv::Mat(static_cast<int>(normImg.height), static_cast<int>(normImg.width), CV_8UC1);
	memcpy(grayImage.data, normImg.buffer.data(), normImg.buffer.size());
	cv::Mat binaryImage = cv::Mat::zeros(static_cast<int>(normImg.height), static_cast<int>(normImg.width), CV_8UC1);

	unsigned int mapWidth = OpenLQM::Core::IntDiv(normImg.width, 4U, true);
	unsigned int mapHeight = OpenLQM::Core::IntDiv(normImg.height, 4U, true);
	unsigned int mapSize = mapWidth * mapHeight;
	cv::Mat dirMap = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);

	if (inputQualityMap) {
		if (qualMap.buffer.size() != mapSize) {
			throw std::invalid_argument(std::string("The input quality map is incompatible with the input image. The quality map must have dimensions ") + std::to_string(mapWidth) + " x " + std::to_string(mapHeight));
		}
	} else {
		if (qualMap.buffer.size() > mapSize) {
			qualMap.buffer.resize(static_cast<std::size_t>(mapSize));
		}
		for (auto& val : qualMap.buffer) {
			val = 0;
		}
		qualMap.buffer.resize(static_cast<std::size_t>(mapSize));
		qualMap.width = mapWidth;
		qualMap.height = mapHeight;
	}
	int maxMinutiae = 1000;
	std::vector<OpenLQM::Core::MinutiaOut> minArray(static_cast<std::size_t>(maxMinutiae));
	cv::Mat gray10pctMap = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);
	cv::Mat gray90pctMap = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);
	cv::Mat grayMedianMap = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);
	cv::Mat grayRangeMap = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);
	cv::Mat grayCountMap = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);
	std::vector<double> maxMagnitude(mapSize);
	std::vector<double> normMagnitude(mapSize);
	std::vector<double> lowFreqMagnitude(mapSize);
	cv::Mat directionChangeMap = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);
	cv::Mat validNeighbors = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);
	cv::Mat curvatureMap = cv::Mat::zeros(static_cast<int>(mapWidth), static_cast<int>(mapHeight), CV_8UC1);

	char* pQualMap = inputQualityMap ? nullptr : reinterpret_cast<char*>(qualMap.buffer.data());
	int numFound = OpenLQM::Core::EmbedMinutiaeDetail(
		nullptr,
		grayImage.data,
		binaryImage.data,
		static_cast<int>(static_cast<unsigned int>(inpImg.resolution)),
		static_cast<int>(normImg.width), static_cast<int>(normImg.height),
		reinterpret_cast<char*>(dirMap.data),
		pQualMap,
		maxMinutiae,
		minArray.data(),
		gray10pctMap.data,
		gray90pctMap.data,
		grayMedianMap.data,
		grayRangeMap.data,
		grayCountMap.data,
		maxMagnitude.data(),
		normMagnitude.data(),
		lowFreqMagnitude.data(),
		directionChangeMap.data,
		validNeighbors.data,
		curvatureMap.data
	);

	if (justQualityMap) {
		return true;
	}

	std::vector<unsigned char> clarityMapNonLca(mapSize);
	std::vector<unsigned char> locq(mapSize);
	OpenLQM::Core::ProcessClarityMap(
		clarityMapNonLca.data(), static_cast<int>(mapWidth), static_cast<int>(mapHeight),
		locq.data(),
		grayRangeMap.data,
		curvatureMap.data,
		directionChangeMap.data,
		lowFreqMagnitude.data(),
		maxMagnitude.data(),
		grayMedianMap.data,
		grayCountMap.data,
		validNeighbors.data,
		normMagnitude.data(),
		qualMap.buffer.data(),
		&inpImg.roi
	);

	if (pFeatureMaps) {
		OpenLQM::Core::CopyMatToFeatureMap(curvatureMap, pFeatureMaps->curvature);
		OpenLQM::Core::CopyMatToFeatureMap(directionChangeMap, pFeatureMaps->direction);
		OpenLQM::Core::CopyMatToFeatureMap(grayCountMap, pFeatureMaps->grayCount);
		OpenLQM::Core::CopyMatToFeatureMap(grayMedianMap, pFeatureMaps->grayMedian);
		OpenLQM::Core::CopyMatToFeatureMap(grayRangeMap, pFeatureMaps->grayRange);
		OpenLQM::Core::CopyVectorToFeatureMap(lowFreqMagnitude, mapWidth, mapHeight, pFeatureMaps->lowFreqMag);
		OpenLQM::Core::CopyVectorToFeatureMap(maxMagnitude, mapWidth, mapHeight, pFeatureMaps->maxMag);
		OpenLQM::Core::CopyVectorToFeatureMap(normMagnitude, mapWidth, mapHeight, pFeatureMaps->normMag);
		OpenLQM::Core::CopyMatToFeatureMap(validNeighbors, pFeatureMaps->validNeighbors);
		pFeatureMaps->quality = qualMap;

		return true;
	}

	return OpenLQM::Core::CalculateMetrics(
		inpImg.buffer.data(),
		static_cast<int>(inpImg.width), static_cast<int>(inpImg.height), static_cast<int>(static_cast<unsigned int>(inpImg.resolution)),
		locq.data(), static_cast<int>(mapWidth), static_cast<int>(mapHeight),
		minArray.data(), numFound,
		metrics,
		pQualMeasures
	);
}

void OPENLQM_API_IMPL OpenLQM::GetAllMetricsFromFingerprint(const OpenLQM::Fingerprint& inpImg, OpenLQM::Metrics& metrics) {
	OpenLQM::QualityMap qualMap;
	GetAllMetricsAndQualityMapFromFingerprint_Base(inpImg, metrics, qualMap, false, false);
}

int OPENLQM_API_IMPL OpenLQM::GetLatentQualityFromFingerprint(const OpenLQM::Fingerprint& inpImg) {
	OpenLQM::Metrics metrics;
	OpenLQM::QualityMap qualMap;
	GetAllMetricsAndQualityMapFromFingerprint_Base(inpImg, metrics, qualMap, false, false);

	return metrics.overallQuality;
}

void OPENLQM_API_IMPL OpenLQM::GetQualityMapFromFingerprint(const OpenLQM::Fingerprint& inpImg, OpenLQM::QualityMap& qualMap) {
	Metrics metrics;
	GetAllMetricsAndQualityMapFromFingerprint_Base(inpImg, metrics, qualMap, false, true);
}

void OPENLQM_API_IMPL OpenLQM::GetAllMetricsAndQualityMapFromFingerprint(const OpenLQM::Fingerprint& inpImg, OpenLQM::Metrics& metrics, OpenLQM::QualityMap& qualMap) {
	GetAllMetricsAndQualityMapFromFingerprint_Base(inpImg, metrics, qualMap, false, false);
}

void OPENLQM_API_IMPL OpenLQM::GetAllMetricsFromFingerprintAndQualityMap(const OpenLQM::Fingerprint& inpImg, OpenLQM::Metrics& metrics, const OpenLQM::QualityMap& qualMap) {
	GetAllMetricsAndQualityMapFromFingerprint_Base(inpImg, metrics, const_cast<OpenLQM::QualityMap&>(qualMap), true, false);
}

void OPENLQM_API_IMPL OpenLQM::ConvertQualityMapCoordinateToFingerprintCoordinate(const OpenLQM::Coordinate& qualCoord, OpenLQM::Coordinate& fingerprintCoord, OpenLQM::PixelDensity fingerprintResolution) {
	int scaleFactor = static_cast<int>(static_cast<unsigned int>(fingerprintResolution) / 125U);
	fingerprintCoord.x = qualCoord.x * scaleFactor;
	fingerprintCoord.y = qualCoord.y * scaleFactor;
}

void OPENLQM_API_IMPL OpenLQM::ConvertFingerprintCoordinateToQualityMapCoordinate(const OpenLQM::Fingerprint& fingerprint, const OpenLQM::Coordinate& fingerprintCoord, OpenLQM::Coordinate& qualCoord) {
	unsigned int ppi = static_cast<unsigned int>(fingerprint.resolution);
	unsigned int resFactor = ppi / 500u;

	unsigned int scaledWidth = std::lrint(static_cast<double>(fingerprint.width) / static_cast<double>(resFactor));
	unsigned int scaledHeight = std::lrint(static_cast<double>(fingerprint.height) / static_cast<double>(resFactor));

	scaledWidth = OpenLQM::Core::WordAlign(scaledWidth, false);
	scaledHeight = OpenLQM::Core::WordAlign(scaledHeight, false);

	scaledWidth = OpenLQM::Core::IntDiv(scaledWidth, 4u, true);
	scaledHeight = OpenLQM::Core::IntDiv(scaledHeight, 4u, true);

	resFactor *= 4u;
	qualCoord.x = std::min(fingerprintCoord.x / static_cast<int>(resFactor), static_cast<int>(scaledWidth - 1u));
	qualCoord.y = std::min(fingerprintCoord.y / static_cast<int>(resFactor), static_cast<int>(scaledHeight - 1u));
}

void OPENLQM_API_IMPL OpenLQM::Supplement::GetAllLowLevelClarityFeatureMapsFromFingerprint(const OpenLQM::Fingerprint& inpImg, OpenLQM::Supplement::FeatureMaps& featureMaps) {
	Metrics metrics;
	QualityMap qualMap;
	GetAllMetricsAndQualityMapFromFingerprint_Base(inpImg, metrics, qualMap, false, false, nullptr, &featureMaps);
}

bool OPENLQM_API_IMPL OpenLQM::Supplement::IsImageNonFingerprint(const OpenLQM::Fingerprint& inpImg) {
	Metrics metrics;
	QualityMap qualMap;
	return !GetAllMetricsAndQualityMapFromFingerprint_Base(inpImg, metrics, qualMap, false, false);
}

void OPENLQM_API_IMPL OpenLQM::Supplement::GetQualityMeasuresFromFingerprintAndQualityMap(const OpenLQM::Fingerprint& inpImg, const OpenLQM::QualityMap& qualMap, OpenLQM::Supplement::QualityMeasures& qualityMeasures) {
	Metrics metrics;
	GetAllMetricsAndQualityMapFromFingerprint_Base(inpImg, metrics, const_cast<OpenLQM::QualityMap&>(qualMap), true, false, &qualityMeasures, nullptr);
}

int OPENLQM_API_IMPL OpenLQM::Supplement::GetLatentQualityFromQualityMeasures(const OpenLQM::Supplement::QualityMeasures& qualityMeasures) {
	if (!OpenLQM::Core::CheckHeuristics(qualityMeasures.LCA1plus, qualityMeasures.LCA3plus, qualityMeasures.Min2plus)) {
		return 0;
	}

	double predictedMatchScore = 7666.16 + (qualityMeasures.Area3plusDiv2 * 32.38) + (qualityMeasures.Area3plus * 17.71) + (qualityMeasures.GFA3plus * 20.87) - (static_cast<double>(qualityMeasures.Min2) * 21.05) + (static_cast<double>(qualityMeasures.Min4) * 7.71);
	int lqMetric12 = static_cast<int>(std::lrint(OpenLQM::Core::CalculateLQMetric12(predictedMatchScore) * 100.0));
	return std::min(lqMetric12, 99);
}

bool OPENLQM_API_IMPL OpenLQM::Supplement::LoadModel(std::string modelPath_) {
	cv::Ptr<cv::ml::RTrees> pNew;
	try {
		pNew = cv::ml::RTrees::load(modelPath_);
	} catch (...) {
		return false;
	}

	const std::vector<int>& roots = pNew->getRoots();
	const std::vector<cv::ml::DTrees::Node>& nodes = pNew->getNodes();
	const std::vector<cv::ml::DTrees::Split>& splits = pNew->getSplits();
	int numClasses = OpenLQM::Impl::GetClassCount(nodes);
	if (numClasses <= 0) {
		throw std::invalid_argument("LoadModel(): Could not detect class count of forest model");
	}


	OpenLQM::Impl::modelPath = modelPath_;
	OpenLQM::Impl::roots = roots;
	OpenLQM::Impl::nodes = nodes;
	OpenLQM::Impl::splits = splits;
	OpenLQM::Impl::numTrees = static_cast<int>(roots.size());
	OpenLQM::Impl::numClasses = numClasses;
	OpenLQM::Impl::loadedModel = pNew;

	return true;
}
