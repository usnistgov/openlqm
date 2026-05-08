/*
 * Copyright 2025 Noblis, Inc.
 * Licensed under the Apache License, Version 2.0.
 *
 * Resolution metadata readers built on libjpeg, libpng, libtiff, and openjpeg —
 * the same families of codecs OpenCV uses for imdecode, without linking FreeImage.
 */

#include "image_resolution.hpp"

#include <csetjmp>

extern "C" {
#include <stdio.h>
#include <jpeglib.h>
#include <png.h>
#include <tiff.h>
#include <tiffio.h>
}

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace {

constexpr float kInchesPerMeter = 39.37007874015748F;

unsigned int DpiToDotsPerMeter(double dpi) {
	if (!(dpi > 0.0) || std::isnan(dpi) || std::isinf(dpi)) {
		return 0U;
	}
	const double dpm = dpi * static_cast<double>(kInchesPerMeter);
	if (dpm <= 0.0 || dpm > static_cast<double>(std::numeric_limits<unsigned int>::max())) {
		return 0U;
	}
	return static_cast<unsigned int>(std::lround(dpm));
}

unsigned int DpcmToDotsPerMeter(double dpcm) {
	if (!(dpcm > 0.0) || std::isnan(dpcm) || std::isinf(dpcm)) {
		return 0U;
	}
	const double dpm = dpcm * 100.0;
	if (dpm <= 0.0 || dpm > static_cast<double>(std::numeric_limits<unsigned int>::max())) {
		return 0U;
	}
	return static_cast<unsigned int>(std::lround(dpm));
}

bool IsJp2Signature(const unsigned char* d, std::size_t n) {
	static const unsigned char soc[12] = {
		0x00, 0x00, 0x00, 0x0C, 0x6A, 0x50, 0x20, 0x20, 0x0D, 0x0A, 0x87, 0x0A};
	return n >= 12 && std::memcmp(d, soc, 12) == 0;
}

bool IsJpegSignature(const unsigned char* d, std::size_t n) {
	return n >= 2 && d[0] == 0xFF && d[1] == 0xD8;
}

bool IsPngSignature(const unsigned char* d, std::size_t n) {
	static const unsigned char sig[8] = {0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A};
	return n >= 8 && std::memcmp(d, sig, 8) == 0;
}

bool IsTiffSignature(const unsigned char* d, std::size_t n) {
	return n >= 4 && ((d[0] == 'I' && d[1] == 'I' && d[2] == 0x2A && d[3] == 0x00)
		|| (d[0] == 'M' && d[1] == 'M' && d[2] == 0x00 && d[3] == 0x2A));
}

bool IsBmpSignature(const unsigned char* d, std::size_t n) {
	return n >= 2 && d[0] == 'B' && d[1] == 'M';
}

std::uint16_t ReadBe16(const unsigned char* p) {
	return static_cast<std::uint16_t>((static_cast<unsigned int>(p[0]) << 8) | static_cast<unsigned int>(p[1]));
}

std::uint32_t ReadLe32(const unsigned char* p) {
	return (static_cast<std::uint32_t>(p[0]))
		| (static_cast<std::uint32_t>(p[1]) << 8)
		| (static_cast<std::uint32_t>(p[2]) << 16)
		| (static_cast<std::uint32_t>(p[3]) << 24);
}

// JP2 capture resolution box ('resc'): ISO/IEC 15444-1 Annex I.7.3.6.1
// Values follow NIST libbiomeval parse_res (PPCM horizontal → dots/m).
unsigned int ReadJp2DotsPerMeterX(const unsigned char* data, std::size_t size) {
	static const unsigned char tag[4] = {'r', 'e', 's', 'c'};
	for (std::size_t i = 0; i + 4 + 10 <= size; ++i) {
		if (std::memcmp(data + i, tag, 4) != 0) {
			continue;
		}
		const unsigned char* p = data + i + 4;
		const std::uint16_t hrN = ReadBe16(p + 4);
		const std::uint16_t hrD = ReadBe16(p + 6);
		const auto hrE = static_cast<int8_t>(p[9]);
		if (hrD == 0) {
			return 0U;
		}
		const double hrPpcm =
			(static_cast<double>(hrN) / static_cast<double>(hrD)) * std::pow(10.0, static_cast<double>(hrE)) / 100.0;
		return DpcmToDotsPerMeter(hrPpcm);
	}
	return 0U;
}

struct JpegErrorMgr {
	jmp_buf setjmp_buffer;
	struct jpeg_error_mgr pub;
};

void JpegErrorExit(j_common_ptr cinfo) {
	auto* err = reinterpret_cast<JpegErrorMgr*>(
		reinterpret_cast<char*>(cinfo->err) - offsetof(JpegErrorMgr, pub));
	longjmp(err->setjmp_buffer, 1);
}

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)  // all locals are trivially destructible; longjmp is safe
#endif
unsigned int ReadJpegDotsPerMeterX(const unsigned char* data, std::size_t size) {
	if (size > static_cast<std::size_t>(std::numeric_limits<unsigned long>::max())) {
		return 0U;
	}
	jpeg_decompress_struct cinfo{};
	JpegErrorMgr jerr{};
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = JpegErrorExit;
	if (setjmp(jerr.setjmp_buffer)) {
		jpeg_destroy_decompress(&cinfo);
		return 0U;
	}
	jpeg_create_decompress(&cinfo);
	jpeg_mem_src(&cinfo, data, static_cast<unsigned long>(size));
	if (jpeg_read_header(&cinfo, TRUE) != JPEG_HEADER_OK) {
		jpeg_destroy_decompress(&cinfo);
		return 0U;
	}
	unsigned int out = 0U;
	if (cinfo.density_unit == 1 && cinfo.X_density > 0) {
		out = DpiToDotsPerMeter(static_cast<double>(cinfo.X_density));
	} else if (cinfo.density_unit == 2 && cinfo.X_density > 0) {
		out = DpcmToDotsPerMeter(static_cast<double>(cinfo.X_density));
	}
	jpeg_destroy_decompress(&cinfo);
	return out;
}
#ifdef _MSC_VER
#pragma warning(pop)
#endif

struct PngMemReader {
	const unsigned char* data = nullptr;
	std::size_t size = 0;
	std::size_t offset = 0;
};

void PngReadCallback(png_structp png, png_bytep out, png_size_t length) {
	auto* st = static_cast<PngMemReader*>(png_get_io_ptr(png));
	if (!st || st->offset > st->size) {
		png_error(png, "read past end");
		return;
	}
	const std::size_t remain = st->size - st->offset;
	if (static_cast<std::size_t>(length) > remain) {
		png_error(png, "read past end");
		return;
	}
	std::memcpy(out, st->data + st->offset, length);
	st->offset += static_cast<std::size_t>(length);
}

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)
#endif
unsigned int ReadPngDotsPerMeterX(const unsigned char* data, std::size_t size) {
	if (size < 8) {
		return 0U;
	}
	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png) {
		return 0U;
	}
	png_infop info = png_create_info_struct(png);
	if (!info) {
		png_destroy_read_struct(&png, nullptr, nullptr);
		return 0U;
	}
	if (setjmp(png_jmpbuf(png))) {
		png_destroy_read_struct(&png, &info, nullptr);
		return 0U;
	}
	PngMemReader mem{data, size, 0};
	png_set_read_fn(png, &mem, PngReadCallback);
	png_read_info(png, info);
	int unit = PNG_RESOLUTION_UNKNOWN;
	png_uint_32 x_pixels_per_meter = 0;
	png_uint_32 y_pixels_per_meter = 0;
	if (png_get_valid(png, info, PNG_INFO_pHYs)) {
		png_get_pHYs(png, info, &x_pixels_per_meter, &y_pixels_per_meter, &unit);
	}
	png_destroy_read_struct(&png, &info, nullptr);
	if (unit != PNG_RESOLUTION_METER || x_pixels_per_meter == 0) {
		return 0U;
	}
	return static_cast<unsigned int>(x_pixels_per_meter);
}
#ifdef _MSC_VER
#pragma warning(pop)
#endif

struct TiffMemState {
	const unsigned char* base = nullptr;
	std::size_t size = 0;
	std::size_t off = 0;
};

static tsize_t TiffRead(thandle_t h, void* buf, tsize_t n) {
	auto* st = static_cast<TiffMemState*>(h);
	if (!st->base || st->off > st->size) {
		return 0;
	}
	const std::size_t remain = st->size - st->off;
	const std::size_t take = static_cast<std::size_t>(n) <= remain ? static_cast<std::size_t>(n) : remain;
	std::memcpy(buf, st->base + st->off, take);
	st->off += take;
	return static_cast<tsize_t>(take);
}

static tsize_t TiffWrite(thandle_t, void*, tsize_t) {
	return 0;
}

static toff_t TiffSeek(thandle_t h, toff_t off, int whence) {
	auto* st = static_cast<TiffMemState*>(h);
	toff_t pos = static_cast<toff_t>(st->off);
	switch (whence) {
	case SEEK_SET:
		pos = off;
		break;
	case SEEK_CUR:
		pos += off;
		break;
	case SEEK_END:
		pos = static_cast<toff_t>(st->size) + off;
		break;
	default:
		return static_cast<toff_t>(-1);
	}
	if (pos < 0 || static_cast<std::uint64_t>(pos) > st->size) {
		return static_cast<toff_t>(-1);
	}
	st->off = static_cast<std::size_t>(pos);
	return pos;
}

static int TiffClose(thandle_t) {
	return 0;
}

static toff_t TiffSize(thandle_t h) {
	return static_cast<toff_t>(static_cast<TiffMemState*>(h)->size);
}

static int TiffMap(thandle_t, void**, toff_t*) {
	return 0;
}

static void TiffUnmap(thandle_t, void*, toff_t) {
}

unsigned int ReadTiffDotsPerMeterX(const unsigned char* data, std::size_t size) {
	TiffMemState st{data, size, 0};
	TIFF* tif = TIFFClientOpen(
		"memory",
		"r",
		reinterpret_cast<thandle_t>(&st),
		TiffRead,
		TiffWrite,
		TiffSeek,
		TiffClose,
		TiffSize,
		TiffMap,
		TiffUnmap);
	if (!tif) {
		return 0U;
	}
	float xres = 0.F;
	uint16_t resunit = RESUNIT_NONE;
	const int hasX = TIFFGetField(tif, TIFFTAG_XRESOLUTION, &xres);
	TIFFGetFieldDefaulted(tif, TIFFTAG_RESOLUTIONUNIT, &resunit);
	TIFFClose(tif);
	if (!hasX || !(xres > 0.F)) {
		return 0U;
	}
	if (resunit == RESUNIT_INCH) {
		return DpiToDotsPerMeter(static_cast<double>(xres));
	}
	if (resunit == RESUNIT_CENTIMETER) {
		return DpcmToDotsPerMeter(static_cast<double>(xres));
	}
	return 0U;
}

unsigned int ReadBmpDotsPerMeterX(const unsigned char* data, std::size_t size) {
	// BITMAPFILEHEADER (14) + DIB header; biXPelsPerMeter only valid for BITMAPINFOHEADER+ (biSize >= 40).
	if (size < 54) {
		return 0U;
	}
	const std::uint32_t dibSize = ReadLe32(data + 14);
	if (dibSize < 40) {
		return 0U;
	}
	const std::int32_t ppm = static_cast<std::int32_t>(ReadLe32(data + 38));
	if (ppm <= 0) {
		return 0U;
	}
	return static_cast<unsigned int>(ppm);
}

} // namespace

namespace OpenLQM::Detail {

unsigned int ReadDotsPerMeterX(const unsigned char* data, std::size_t size) {
	if (!data || size < 12) {
		return 0U;
	}
	if (IsJp2Signature(data, size)) {
		return ReadJp2DotsPerMeterX(data, size);
	}
	if (IsJpegSignature(data, size)) {
		return ReadJpegDotsPerMeterX(data, size);
	}
	if (IsPngSignature(data, size)) {
		return ReadPngDotsPerMeterX(data, size);
	}
	if (IsTiffSignature(data, size)) {
		return ReadTiffDotsPerMeterX(data, size);
	}
	if (IsBmpSignature(data, size)) {
		return ReadBmpDotsPerMeterX(data, size);
	}
	return 0U;
}

} // namespace OpenLQM::Detail
