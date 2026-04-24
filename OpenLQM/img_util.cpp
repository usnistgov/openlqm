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

#include <openlqm/openlqm_img_util.hpp>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <cstring>

namespace OpenLQM {
    namespace Core {
        int GetMaxPadding(const int map_windowsize, const int map_windowoffset, const int dirbin_grid_w, const int dirbin_grid_h) {
            int dft_pad, dirbin_pad, max_pad;
            double diag;
            double pad;


            /* 1. Compute pad required for rotated windows used in DFT analyses. */

            /* Explanation of DFT padding:

                            B---------------------
                            |      window        |
                            |                    |
                            |                    |
                            |      A.......______|__________
                            |      :      :      |
                            |<-C-->: block:      |
                        <--|--D-->:      :      | image
                            |      ........      |
                            |      |             |
                            |      |             |
                            |      |             |
                            ----------------------
                                    |
                                    |
                                    |

                Pixel A = Origin of entire fingerprint image
                        = Also origin of first block in image. Each pixel in
                            this block gets the same DFT results computed from
                            the surrounding window.  Note that in general
                            blocks are adjacent and non-overlapping.

                Pixel B = Origin of surrounding window in which DFT
                            analysis is conducted.  Note that this window is not
                            completely contained in the image but extends to the
                            top and to the right.

                Distance C = Number of pixels in which the window extends
                            beyond the image (map_windowoffset).

                Distance D = Amount of padding required to hold the entire
                            rotated window in memory.

            */

            /* Compute pad as difference between the MAP windowsize           */
            /* and the diagonal distance of the window.                       */
            /* (DFT grids are computed with pixel offsets RELATIVE2ORIGIN.)   */
            diag = std::sqrt(static_cast<double>(2.0 * map_windowsize * map_windowsize));
            pad = (diag-map_windowsize)/static_cast<double>(2.0);
            /* Need to truncate precision so that answers are consistent  */
            /* on different computer architectures when rounding doubles. */
            pad = trunc_dbl_precision(pad);
            /* Must add the window offset to the rotational padding. */
            dft_pad = sround(pad) + map_windowoffset;

            /* 2. Compute pad required for rotated blocks used in directional  */
            /*    binarization.  Binarization blocks are applied to each pixel */
            /*    in the input image.                                          */
            diag = sqrt(static_cast<double>((dirbin_grid_w*dirbin_grid_w)+
                                (dirbin_grid_h*dirbin_grid_h)));
            /* Assumption: all grid centers reside in valid/allocated memory. */
            /* (Dirbin grids are computed with pixel offsets RELATIVE2CENTER.) */
            pad = (diag-1)/static_cast<double>(2.0);
            /* Need to truncate precision so that answers are consistent */
            /* on different computer architectures when rounding doubles. */
            pad = trunc_dbl_precision(pad);
            dirbin_pad = sround(pad);

            max_pad = std::max(dft_pad, dirbin_pad);

            /* Return the maximum of the two required paddings.  This padding will */
            /* be sufficiently large for all purposes, so that padding of the      */
            /* input image will only be required once.                             */
            return max_pad;
        }

        void PadUCharImage(std::vector<unsigned char>& out_paddedData, int& ow, int& oh,
                            unsigned char *idata, const int iw, const int ih,
                            const int pad, const int pad_value)
        {
            unsigned char *pptr;
            int pw, ph;
            int pad2, psize;

            /* Account for pad on both sides of image */
            pad2 = pad<<1;

            /* Compute new pad sizes */
            pw = iw + pad2;
            ph = ih + pad2;
            psize = pw * ph;

            /* Allocate padded image */
            out_paddedData.resize(static_cast<std::size_t>(psize));
            unsigned char* pdata = out_paddedData.data();

            /* Initialize values to a constant PAD value */
            memset(pdata, pad_value, static_cast<std::size_t>(psize));

            /* Copy input image into padded image one scanline at a time */
            pptr = pdata + (pad * pw) + pad;
            for (int i = 0; i < ih; i++) {
                memcpy(pptr, idata, static_cast<std::size_t>(iw));
                idata += iw;
                pptr += pw;
            }

            ow = pw;
            oh = ph;
        }

        void Bits8to6(unsigned char *idata, const int iw, const int ih) {
            int isize = iw * ih;
            for (int i = 0; i < isize; i++) {
                /* Divide every pixel value by 4 so that [0..256) -> [0..64) */
                *idata++ >>= 2;
            }
        }

        void GenImageMaps(std::vector<int>& out_dmap, std::vector<int>& out_lcmap, std::vector<int>& out_lfmap, std::vector<int>& out_hcmap,
                        int *omw, int *omh,
                        unsigned char *pdata, const int pw, const int ph,
                        const Dir2Rad& dir2rad, const DFTWaves& dftwaves,
                        const RotGrids& dftGrids, const LQMParams& lqmParams,
                        unsigned char *gray10pctMap, unsigned char *gray90pctMap, unsigned char *grayMedianMap, unsigned char *grayRangeMap, unsigned char *grayCountMap,
                        double *maxMagnitude, double *normMagnitude, double *lowFreqMagnitude,
                        unsigned char *directionChangeMap,
                        unsigned char *validNeighbors, unsigned char *curvatureMap) {

            /* Compute block offsets for the entire image, accounting for pad */
            /* BlockOffsets() assumes square block (grid)                     */
            if (dftGrids.grid_w != dftGrids.grid_h) {
                throw std::invalid_argument("GenImageMaps: DFT grids must be square");
            }

            /* Compute unpadded image dimensions. */
            std::vector<int> blockOffsets;
            int iw = pw - (dftGrids.pad<<1);
            int ih = ph - (dftGrids.pad<<1);
            int mw = 0, mh = 0;
            BlockOffsets(blockOffsets, mw, mh, iw, ih,
                        dftGrids.pad, lqmParams.blocksize);

            /* Generate initial Direction Map and Low Contrast Map */
            GenInitialMaps(out_dmap, out_lcmap,
                            out_lfmap, blockOffsets.data(), mw, mh,
                            pdata, pw, ph, dftwaves, dftGrids, lqmParams,
                            gray10pctMap, gray90pctMap, grayMedianMap, grayRangeMap, grayCountMap,
                            maxMagnitude, normMagnitude, lowFreqMagnitude);

            // Retain initial Direction Map
            std::vector<int> initialDirectionMap(out_dmap);

            /* Blur low flow map */
            MorphTFMap(out_lfmap.data(), mw, mh);

            /* Remove directions that are inconsistent with neighbors */
            RemoveInconDirs(out_dmap.data(), mw, mh, dir2rad, lqmParams);

            /* Smooth Direction Map values with their neighbors */
            SmoothDirectionMap(out_dmap.data(), out_lcmap.data(), mw, mh,
                                    dir2rad, lqmParams);

            /* Interpolate INVALID direction blocks with their valid neighbors. */
            InterpolateDirectionMap(out_dmap.data(), out_lcmap.data(),
                                                mw, mh, lqmParams);

            /* May be able to skip steps 6 and/or 7 if computation time */
            /* is a critical factor.                                    */

            /* Remove directions that are inconsistent with neighbors */
            RemoveInconDirs(out_dmap.data(), mw, mh, dir2rad, lqmParams);

            /* Smooth Direction Map values with their neighbors. */
            SmoothDirectionMap(out_dmap.data(), out_lcmap.data(), mw, mh,
                                    dir2rad, lqmParams);

            /* Set the Direction Map values in the image margin to INVALID. */
            SetMarginBlocks(out_dmap.data(), mw, mh, LQM_INVALID_DIR);

            /* Generate High Curvature Map from interpolated Direction Map. */
            GenHighCurveMap(out_hcmap, out_dmap.data(), mw, mh,
                                            lqmParams, validNeighbors, curvatureMap);

            GenDirectionChangeMap(directionChangeMap, initialDirectionMap.data(), out_dmap.data(), mw, mh);

            *omw = mw;
            *omh = mh;
        }

        void BlockOffsets(std::vector<int>& out_blockOffsets, int& ow, int& oh,
                const int iw, const int ih, const int pad, const int blocksize) {

            /* Test if unpadded image is smaller than a single block */
            if ((iw < blocksize) || (ih < blocksize)) {
                std::string blocksizeStr = std::to_string(blocksize);
                throw std::invalid_argument(std::string("BlockOffsets: Image must be at least ") + blocksizeStr + " by " + blocksizeStr + " in size");
            }

            /* Compute padded width and height of image */
            int pad2 = pad<<1;
            int pw = iw + pad2;

            /* Compute the number of columns and rows of blocks in the image. */
            /* Take the ceiling to account for "leftovers" at the right and   */
            /* bottom of the unpadded image */
            int bw = static_cast<int>(std::ceil(iw / static_cast<double>(blocksize)));
            int bh = static_cast<int>(std::ceil(ih / static_cast<double>(blocksize)));

            /* Total number of blocks in the image */
            int bsize = bw*bh;

            /* The index of the last column */
            int lastbw = bw - 1;
            /* The index of the last row */
            int lastbh = bh - 1;

            /* Allocate list of block offsets */
            out_blockOffsets.resize(static_cast<std::size_t>(bsize));

            /* Current block index */
            int bi = 0;

            /* Current offset from top of padded image to start of new row of  */
            /* unpadded image blocks. It is initialize to account for the      */
            /* padding and will always be indented the size of the padding     */
            /* from the left edge of the padded image.                         */
            int blkrow_start = (pad * pw) + pad;

            /* Number of pixels in a row of blocks in the padded image */
            int blkrow_size = pw * blocksize;  /* row width X block height */

            /* Foreach non-overlapping row of blocks in the image */
            for (int by = 0; by < lastbh; by++, blkrow_start += blkrow_size) {
                /* `offset` is the current offset from the top of the padded image to beginning of */
                /* the next block */
                int offset = blkrow_start;

                /* Foreach non-overlapping column of blocks in the image */
                for (int bx = 0; bx < lastbw; bx++) {
                    /* Store current block offset */
                    out_blockOffsets[static_cast<std::size_t>(bi++)] = offset;
                    /* Bump to the beginning of the next block */
                    offset += blocksize;
                }

                /* Compute and store "left-over" block in row.    */
                /* This is the block in the last column of row.   */
                /* Start at far right edge of unpadded image data */
                /* and come in BLOCKSIZE pixels.                  */
                out_blockOffsets[static_cast<std::size_t>(bi++)] = blkrow_start + iw - blocksize;

                /* Bump to beginning of next row of blocks (performed by the for loop update `blkrow_start += blkrow_size`) */
            }

            /* Compute and store "left-over" row of blocks at bottom of image */
            /* Start at bottom edge of unpadded image data and come up        */
            /* BLOCKSIZE pixels. This too must account for padding.           */
            blkrow_start = ((pad + ih - blocksize) * pw) + pad;

            /* Foreach non-overlapping column of blocks in last row of the image */
            for (int bx = 0, offset = blkrow_start; bx < lastbw; bx++, offset += blocksize) {
                /* Store current block offset */
                out_blockOffsets[static_cast<std::size_t>(bi++)] = offset;
            }

            /* Compute and store last "left-over" block in last row.      */
            /* Start at right edge of unpadded image data and come in     */
            /* BLOCKSIZE pixels.                                          */
            out_blockOffsets[static_cast<std::size_t>(bi++)] = blkrow_start + iw - blocksize;

            ow = bw;
            oh = bh;
        }

        void GenInitialMaps(std::vector<int>& out_dmap, std::vector<int>& out_lcmap, std::vector<int>& out_lfmap,
                        int *blkoffs, const int mw, const int mh,
                        unsigned char *pdata, const int pw, const int ph,
                        const DFTWaves& dftWaves, const RotGrids& dftGrids,
                        const LQMParams& lqmParams,
                        unsigned char *gray10pctMap, unsigned char *gray90pctMap, unsigned char *grayMedianMap, unsigned char *grayRangeMap, unsigned char *grayCountMap,
                        double *maxMagnitude, double *normMagnitude, double *lowFreqMagnitude) {

            /* Compute total number of blocks in map */
            int bsize = mw * mh;

            /* Initialize Direction Map to INVALID (-1) */
            out_dmap = std::vector<int>(static_cast<std::size_t>(bsize), LQM_INVALID_DIR);
            /* Initialize Low Contrast Map to FALSE (0) */
            out_lcmap = std::vector<int>(static_cast<std::size_t>(bsize), 0);
            /* Initialize Low Flow Map to FALSE (0). */
            out_lfmap = std::vector<int>(static_cast<std::size_t>(bsize), 0);

            /* Allocate DFT directional power vectors */
            std::vector<std::vector<double>> powers;
            AllocDirPowers(powers, dftWaves.nwaves, dftGrids.ngrids);

            /* Allocate DFT power statistic arrays */
            /* Compute length of statistics arrays.  Statistics not needed   */
            /* for the first DFT wave, so the length is number of waves - 1. */
            int nstats = dftWaves.nwaves - 1;
            std::vector<int> wis;
            std::vector<double> powmaxs;
            std::vector<int> powmaxDirs;
            std::vector<double> pownorms;
            AllocPowerStats(wis, powmaxs, powmaxDirs, pownorms, nstats);

            /* Compute special window origin limits for determining low contrast.  */
            /* These pixel limits avoid analyzing the padded borders of the image. */
            int xminlimit = dftGrids.pad;
            int yminlimit = dftGrids.pad;
            int xmaxlimit = pw - dftGrids.pad - lqmParams.windowsize - 1;
            int ymaxlimit = ph - dftGrids.pad - lqmParams.windowsize - 1;

            /* Foreach block in image ... */
            for (int bi = 0; bi < bsize; ++bi) {
                /* Adjust block offset from pointing to block origin to pointing */
                /* to surrounding window origin.                                 */
                int dft_offset = blkoffs[bi] - (lqmParams.windowoffset * pw) - lqmParams.windowoffset;

                /* Compute pixel coords of window origin. */
                int win_x = dft_offset % pw;
                int win_y = static_cast<int>(dft_offset / pw);

                /* Make sure the current window does not access padded image pixels */
                /* for analyzing low contrast.                                      */
                win_x = std::max(xminlimit, win_x);
                win_x = std::min(xmaxlimit, win_x);
                win_y = std::max(yminlimit, win_y);
                win_y = std::min(ymaxlimit, win_y);
                int low_contrast_offset = (win_y * pw) + win_x;

                /* If block is low contrast ... */
                // Sets low_contrast_map as Boolean
                if (LowContrastBlock(low_contrast_offset, lqmParams.windowsize,
                                            pdata, pw, lqmParams,
                                            &(gray10pctMap[bi]), &(gray90pctMap[bi]), &(grayMedianMap[bi]), &(grayRangeMap[bi]), &(grayCountMap[bi]))) {

                    /* Block is low contrast ... */
                    out_lcmap[static_cast<std::size_t>(bi)] = LQM_TRUE;
                    /* Direction Map's block is already set to INVALID. */
                }
                /* Otherwise, sufficient contrast for DFT processing ... */
                else {
                    /* Compute DFT powers */
                    DFTDirPowers(powers, pdata, low_contrast_offset,
                                        dftWaves, dftGrids);

                    /* Compute DFT power statistics, skipping first applied DFT  */
                    /* wave.  This is dependent on how the primary and secondary */
                    /* direction tests work below.                               */
                    DFTPowerStats(wis.data(), powmaxs.data(), powmaxDirs.data(), pownorms.data(), powers,
                                            1, dftWaves.nwaves, dftGrids.ngrids);

                    /* Conduct primary direction test and return powmaxs, pownorms */
                    int blkdir = PrimaryDirTest(powers, wis.data(), powmaxs.data(), powmaxDirs.data(),
                                            pownorms.data(), nstats, lqmParams,
                                            &(maxMagnitude[static_cast<std::size_t>(bi)]), &(normMagnitude[static_cast<std::size_t>(bi)]), &(lowFreqMagnitude[static_cast<std::size_t>(bi)]));

                    if (blkdir != LQM_INVALID_DIR) {
                        out_dmap[static_cast<std::size_t>(bi)] = blkdir;
                    }
                    else {
                        /* Conduct secondary (fork) direction test */
                        blkdir = SecondaryForkTest(powers, wis.data(), powmaxs.data(), powmaxDirs.data(),
                                            pownorms.data(), lqmParams,
                                            &(maxMagnitude[static_cast<std::size_t>(bi)]), &(normMagnitude[static_cast<std::size_t>(bi)]), &(lowFreqMagnitude[static_cast<std::size_t>(bi)]));
                        if(blkdir != LQM_INVALID_DIR) {
                            out_dmap[static_cast<std::size_t>(bi)] = blkdir;
                        }
                        else {
                            /* The current direction in Direction Map remains LQM_INVALID_DIR */
                            /* Flag the block as having LOW RIDGE FLOW. */
                            out_lfmap[static_cast<std::size_t>(bi)] = LQM_TRUE;
                        }
                    }

                } /* End DFT */
            } /* bi */
        }

        void AllocDirPowers(std::vector<std::vector<double>>& out_powers, const int nwaves, const int ndirs)
        {
            /* Allocate list of lists containing power vectors */
            out_powers.resize(static_cast<std::size_t>(nwaves));

            /* Foreach DFT wave ... */
            for (int w = 0; w < nwaves; w++) {
                /* Allocate power vector for all directions */
                out_powers[static_cast<std::size_t>(w)].resize(static_cast<std::size_t>(ndirs));
            }
        }

        void AllocPowerStats(std::vector<int>& out_wis, std::vector<double>& out_powmaxs, std::vector<int>& out_powmaxDirs,
                            std::vector<double>& out_pownorms, const int nstats) {

            /* Allocate DFT wave index vector */
            out_wis.resize(static_cast<std::size_t>(nstats));

            /* Allocate max power vector */
            out_powmaxs.resize(static_cast<std::size_t>(nstats));

            /* Allocate max power direction vector */
            out_powmaxDirs.resize(static_cast<std::size_t>(nstats));

            /* Allocate normalized power vector */
            out_pownorms.resize(static_cast<std::size_t>(nstats));
        }

        bool LowContrastBlock(const int blkoffset, const int blocksize,
                            unsigned char *pdata, const int pw,
                            const LQMParams& lqmParams,
                            unsigned char *pct10, unsigned char *pct90, unsigned char *median, unsigned char *range, unsigned char *numgrays) {

            int pixtable[LQM_MAX_6_BIT_GRAY] = {};
            int numpix = blocksize * blocksize;

            double tdbl = (lqmParams.percentile_min_max/100.0) * static_cast<double>(numpix-1);
            tdbl = trunc_dbl_precision(tdbl);
            int prctthresh = sround(tdbl);

            unsigned char* pRow = pdata + blkoffset;
            for (int py = 0; py < blocksize; ++py) {
                unsigned char* pCol = pRow;
                for (int px = 0; px < blocksize; ++px) {
                    ++pixtable[*pCol];
                    ++pCol;
                }
                pRow += pw;
            }

            // pixtable now contains a histogram of the 6-bit (0-63) grayscale values

            int pixsum = 0;
            int prctmin = -1;
            int prctmax = -1;
            int medianpix = -1;
            int graycount = 0;

            for (int pi = 0; pi < LQM_MAX_6_BIT_GRAY; ++pi) {
                pixsum += pixtable[pi];
                if (pixtable[pi]) {
                    graycount++;
                }
                if (pixsum >= prctthresh) {
                    if (prctmin<0) {
                        prctmin = pi;
                    } // 10th percentile
                    if (pixsum >= (numpix / 2)) {
                        if (medianpix<0) {
                            medianpix = pi;
                        }
                        if (pixsum >= (numpix - 1 - prctthresh)) {
                            if (prctmax<0) {
                                prctmax = pi;
                            } // 90th percentile
                        }
                    }
                }
            }

            if (prctmin < 0) {
                throw std::invalid_argument("LowContrastBlock: 10th percentile pixel not found");
            }
            if (medianpix < 0) {
                throw std::invalid_argument("LowContrastBlock: Median pixel not found");
            }
            if (prctmax < 0)   {
                throw std::invalid_argument("LowContrastBlock: 90th percentile pixel not found");
            }

            *pct10 =  static_cast<unsigned char> (prctmin*4); // *4 is to adjust 6-bit to 8-bit
            *pct90 =  static_cast<unsigned char> (prctmax*4); // *4 is to adjust 6-bit to 8-bit

            *median = static_cast<unsigned char> (medianpix*4); // adjust 6-bit to 8-bit

            *range = *pct90 - *pct10;
            *numgrays = static_cast<unsigned char> (graycount*4); // *4 is to rescale to 8-bit

            int delta = prctmax - prctmin;

            return delta < lqmParams.min_contrast_delta;
        }

        void DFTDirPowers(std::vector<std::vector<double>>& powers, unsigned char *pdata,
                    const int blkoffset,
                    const DFTWaves& dftWaves, const RotGrids& dftGrids) {

            /* Allocate line sum vector, and initialize to zeros */
            /* This routine requires square block (grid), so ERROR otherwise. */
            if (dftGrids.grid_w != dftGrids.grid_h) {
                throw std::invalid_argument("DFTDirPowers: DFT grids must be square");
            }

            std::vector<int> rowsums(static_cast<std::size_t>(dftGrids.grid_w));

            /* Foreach direction ... */
            for (int dir = 0; dir < dftGrids.ngrids; dir++) {
                /* Compute vector of line sums from rotated grid */
                unsigned char* blkptr = pdata + blkoffset;
                SumRotBlockRows(rowsums, blkptr,
                                    dftGrids.grids[static_cast<std::size_t>(dir)].data(), dftGrids.grid_w);

                /* Foreach DFT wave ... */
                for (int w = 0; w < dftWaves.nwaves; w++) {
                    DFTPower(&(powers[static_cast<std::size_t>(w)][static_cast<std::size_t>(dir)]), rowsums.data(),
                            dftWaves.waves[static_cast<std::size_t>(w)], dftWaves.wavelen);
                }
            }
        }

        void SumRotBlockRows(std::vector<int>& rowsums, const unsigned char *blkptr,
                                const int *grid_offsets, const int blocksize) {

            /* Initialize rotation offset index. */
            int gi = 0;

            /* For each row in block ... */
            for (int iy = 0; iy < blocksize; iy++) {
                /* The sums are accumlated along the rotated rows of the grid, */
                /* so initialize row sum to 0.                                 */
                rowsums[static_cast<std::size_t>(iy)] = 0;
                /* Foreach column in block ... */
                for (int ix = 0; ix < blocksize; ix++) {
                    /* Accumulate pixel value at rotated grid position in image */
                    rowsums[static_cast<std::size_t>(iy)] += *(blkptr + grid_offsets[static_cast<std::size_t>(gi)]);
                    ++gi;
                }
            }
        }

        void DFTPower(double *power, const int *rowsums,
                    const DFTWave& wave, const int wavelen) {

        /* Initialize accumulators */
            double cospart = 0.0;
            double sinpart = 0.0;

            /* Accumulate cos and sin components of DFT. */
            for (int i = 0; i < wavelen; ++i) {
                /* Multiply each rotated row sum by its        */
                /* corresponding cos or sin point in DFT wave. */
                cospart += (rowsums[static_cast<std::size_t>(i)] * wave.cosLut[static_cast<std::size_t>(i)]);
                sinpart += (rowsums[static_cast<std::size_t>(i)] * wave.sinLut[static_cast<std::size_t>(i)]);
            }

            /* Power is the sum of the squared cos and sin components */
            *power = (cospart * cospart) + (sinpart * sinpart);
        }

        void DFTPowerStats(int *wis, double *powmaxs, int *powmax_dirs,
                            double *pownorms, std::vector<std::vector<double>>& powers,
                            const int fw, const int tw, const int ndirs) {

            for (int w = fw, i = 0; w < tw; ++w, ++i) {
                GetMaxNorm(&(powmaxs[static_cast<std::size_t>(i)]), &(powmax_dirs[static_cast<std::size_t>(i)]),
                                &(pownorms[static_cast<std::size_t>(i)]), powers[static_cast<std::size_t>(w)].data(), ndirs);
            }

            /* Get sorted order of applied DFT waves based on normalized power */
            SortDFTWaves(wis, powmaxs, pownorms, tw-fw);
        }

        void GetMaxNorm(double *powmax, int *powmax_dir,
                    double *pownorm, const double *power_vector, const int ndirs) {

            /* Find max power value and store corresponding direction */
            double max_v = power_vector[0];
            int max_i = 0;

            /* Sum the total power in a block at a given direction */
            double powsum = power_vector[0];

            /* For each direction ... */
            for (int dir = 1; dir < ndirs; dir++) {
                powsum += power_vector[dir];
                if (power_vector[dir] > max_v) {
                    max_v = power_vector[dir];
                    max_i = dir;
                }
            }

            *powmax = max_v;
            *powmax_dir = max_i;

            /* Powmean is used as denominator for pownorm, so setting  */
            /* a non-zero minimum avoids possible division by zero.    */
            double powmean = std::max(powsum, LQM_MIN_POWER_SUM)/static_cast<double>(ndirs);

            *pownorm = *powmax / powmean;
        }

        void SortDFTWaves(int *wis, const double *powmaxs, const double *pownorms,
                        const int nstats) {

            std::vector<double> pownorms2(static_cast<std::size_t>(nstats));

            for (int i = 0; i < nstats; ++i) {
                /* Wis will hold the sorted statistic indices when all is done. */
                wis[static_cast<std::size_t>(i)] = i;
                /* This is normalized squared max power. */
                pownorms2[static_cast<std::size_t>(i)] = powmaxs[static_cast<std::size_t>(i)] * pownorms[static_cast<std::size_t>(i)];
            }

            // todo: consider replacing this with a better sort
            /* Sort the statistic indices on the normalized squared power. */
            BubbleSortDoubleDec(pownorms2.data(), wis, nstats);
        }

        void BubbleSortDoubleDec(double *ranks, int *items,  const int len) {
            bool done = false;
            int n = len;
            while (!done) {
                done = true;
                for (int i = 1, p = 0; i < n; ++i, ++p) {
                    /* If previous rank is < current rank ... */
                    if (ranks[p] < ranks[i]) {
                        /* Swap ranks */
                        double trank = ranks[i];
                        ranks[i] = ranks[p];
                        ranks[p] = trank;
                        /* Swap corresponding items */
                        int titem = items[i];
                        items[i] = items[p];
                        items[p] = titem;
                        done = false;
                    }
                }
                --n;
            }
        }

        int PrimaryDirTest(std::vector<std::vector<double>>& powers, const int *wis,
                    const double *powmaxs, const int *powmax_dirs,
                    const double *pownorms, const int nstats,
                    const LQMParams& lqmParams,
                    double *maxMagnitude, double *normMagnitude, double *lowFreqMagnitude) {

            /* Look at max power statistics in decreasing order ... */
            for (int w = 0; w < nstats; ++w) {
                /* 1. Test magnitude of current max power (Ex. Thresh==100000)   */
                if (
                    /* 1. Test magnitude of current max power (Ex. Thresh==100000)   */
                    (powmaxs[static_cast<std::size_t>(wis[static_cast<std::size_t>(w)])] > lqmParams.powmax_min) &&
                    /* 2. Test magnitude of normalized max power (Ex. Thresh==3.8)   */
                    (pownorms[static_cast<std::size_t>(wis[static_cast<std::size_t>(w)])] > lqmParams.pownorm_min) &&
                    /* 3. Test magnitude of power of lowest DFT frequency at current */
                    /* max power direction and make sure it is not too big.          */
                    /* (Ex. Thresh==50000000)                                        */
                    (powers[0][static_cast<std::size_t>(powmax_dirs[static_cast<std::size_t>(wis[static_cast<std::size_t>(w)])])] <= lqmParams.powmax_max)
                ) {

                        /* If ALL 3 criteria met, return current max power direction. */
                        *maxMagnitude = powmaxs[wis[static_cast<std::size_t>(w)]];
                        *normMagnitude = pownorms[wis[static_cast<std::size_t>(w)]];
                        *lowFreqMagnitude = powers[0][static_cast<std::size_t>(powmax_dirs[static_cast<std::size_t>(wis[static_cast<std::size_t>(w)])])];

                        return powmax_dirs[static_cast<std::size_t>(wis[static_cast<std::size_t>(w)])];

                }
            }

            /* Otherwise test failed. */

            *maxMagnitude = 0;
            *normMagnitude = 0;
            *lowFreqMagnitude = 0;

            return LQM_INVALID_DIR;
        }

        int SecondaryForkTest(std::vector<std::vector<double>>& powers, const int *wis,
            const double *powmaxs, const int *powmax_dirs,
            const double *pownorms,
            const LQMParams& lqmParams,
            double *maxMagnitude, double *normMagnitude, double *lowFreqMagnitude)
        {
            /* Relax the normalized power threshold under fork conditions. */
            double fork_pownorm_min = lqmParams.fork_pct_pownorm * lqmParams.pownorm_min;

                /* 1. Test magnitude of largest max power (Ex. Thresh==100000)   */
            if ((powmaxs[wis[0]] > lqmParams.powmax_min) &&
                /* 2. Test magnitude of corresponding normalized power           */
                /*    (Ex. Thresh==2.85)                                         */
                (pownorms[wis[0]] >= fork_pownorm_min) &&
                /* 3. Test magnitude of power of lowest DFT frequency at largest */
                /* max power direction and make sure it is not too big.          */
                /* (Ex. Thresh==50000000)                                        */
                (powers[0][static_cast<std::size_t>(powmax_dirs[wis[0]])] <= lqmParams.powmax_max))
            {
                /* Add FORK_INTERVALs to current direction modulo NDIRS */
                int rdir = (powmax_dirs[wis[0]] + lqmParams.fork_interval) % lqmParams.num_directions;

                /* Subtract FORK_INTERVALs from direction modulo NDIRS  */
                /* For example, FORK_INTERVAL==2 & NDIRS==16, then      */
                /*            ldir = (dir - (16-2)) % 16                */
                /* which keeps result in proper modulo range.           */
                int ldir = (powmax_dirs[wis[0]] + lqmParams.num_directions - lqmParams.fork_interval) % lqmParams.num_directions;

                /* Set forked angle threshold to be a % of the max directional */
                /* power. (Ex. thresh==0.7*powmax)                             */
                double fork_pow_thresh = powmaxs[wis[0]] * lqmParams.fork_pct_powmax;

                /* Look up and test the computed power for the left and right    */
                /* fork directions.s                                             */
                /* The power stats (and thus wis) are on the range [0..nwaves-1) */
                /* as the statistics for the first DFT wave are not included.    */
                /* The original power vectors exist for ALL DFT waves, therefore */
                /* wis indices must be added by 1 before addressing the original */
                /* powers vector.                                                */
                /* OpenLQM permits one and only one of the fork angles to exceed */
                /* the relative power threshold.                                 */
                if (((powers[static_cast<std::size_t>(wis[0]+1)][static_cast<std::size_t>(ldir)] <= fork_pow_thresh) ||
                    (powers[static_cast<std::size_t>(wis[0]+1)][static_cast<std::size_t>(rdir)] <= fork_pow_thresh)) &&
                    ((powers[static_cast<std::size_t>(wis[0]+1)][static_cast<std::size_t>(ldir)] > fork_pow_thresh) ||
                    (powers[static_cast<std::size_t>(wis[0]+1)][static_cast<std::size_t>(rdir)] > fork_pow_thresh)))
                {
                    /* If ALL the above criteria hold, then return the direction */
                    /* of the largest max power.                                 */

                    // Returning this sets only a relative handful of pixels
                    *maxMagnitude = powmaxs[wis[0]];
                    *normMagnitude = pownorms[wis[0]];
                    *lowFreqMagnitude = powers[0][static_cast<std::size_t>(powmax_dirs[wis[0]])];

                    return(powmax_dirs[static_cast<std::size_t>(wis[0])]);
                }
            }

            /* Otherwise test failed. */

            *maxMagnitude = 0;
            *normMagnitude = 0;
            *lowFreqMagnitude = 0;

            return LQM_INVALID_DIR;
        }

        void MorphTFMap(int *tfmap, const int mw, const int mh) {

            /* Convert TRUE/FALSE map into a binary byte image. */
            std::vector<unsigned char> cimage(static_cast<std::size_t>(mw * mh));

            unsigned char* cptr = cimage.data();
            int* mptr = tfmap;
            for (int i = 0; i < mw * mh; ++i) {
                *cptr++ = static_cast<unsigned char>(*mptr++);
            }
            std::vector<unsigned char> mimage(static_cast<std::size_t>(mw * mh));

            DilateUCharImage(cimage.data(), mimage.data(), mw, mh);
            DilateUCharImage(mimage.data(), cimage.data(), mw, mh);
            ErodeUCharImage(cimage.data(), mimage.data(), mw, mh);
            ErodeUCharImage(mimage.data(), cimage.data(), mw, mh);

            cptr = cimage.data();
            mptr = tfmap;
            for (int i = 0; i < mw * mh; ++i) {
                *mptr++ = int(*cptr++);
            }
        }

        void DilateUCharImage(unsigned char *inp, unsigned char *out, const int iw, const int ih) {
            memcpy(out, inp, static_cast<std::size_t>(iw * ih));

            /* for all pixels. set pixel if there is at least one true neighbor */
            unsigned char* itr = inp;
            unsigned char* otr = out;
            for (int row = 0; row < ih; ++row ) {
                for (int col = 0; col < iw; ++col) {
                    if (!*itr)     /* pixel is already true, neighbors irrelevant */
                    {
                        /* more efficient with C's left to right evaluation of     */
                        /* conjuctions. E N S functions not executed if W is false */
                        if (
                            GetWestUChar(itr, col, static_cast<unsigned char>(0)) ||
                            GetEastUChar(itr, col, iw, static_cast<unsigned char>(0)) ||
                            GetNorthUChar(itr, row, iw, static_cast<unsigned char>(0)) ||
                            GetSouthUChar(itr, row, iw, ih, static_cast<unsigned char>(0))
                        ) {
                            *otr = 1;
                        }
                    }
                    itr++ ; otr++;
                }
            }
        }

        void ErodeUCharImage(unsigned char *inp, unsigned char *out, const int iw, const int ih) {
            memcpy(out, inp, static_cast<std::size_t>(iw * ih));

            /* for true pixels. kill pixel if there is at least one false neighbor */
            unsigned char* itr = inp;
            unsigned char* otr = out;
            for (int row = 0; row < ih; ++row ) {
                for (int col = 0; col < iw; ++col) {
                    if (*itr)      /* erode only operates on true pixels */
                    {
                        /* more efficient with C's left to right evaluation of     */
                        /* conjuctions. E N S functions not executed if W is false */
                        if (!(
                            GetWestUChar(itr, col, static_cast<unsigned char>(1)) &&
                            GetEastUChar(itr, col, iw, static_cast<unsigned char>(1)) &&
                            GetNorthUChar(itr, row, iw, static_cast<unsigned char>(1)) &&
                            GetSouthUChar(itr, row, iw, ih, static_cast<unsigned char>(1))
                        )) {
                            *otr = 0;
                        }
                    }
                    itr++ ; otr++;
                }
            }
        }

        unsigned char GetSouthUChar(unsigned char *ptr, const int row, const int iw, const int ih, unsigned char defaultVal)
        {
            /* catch case where image is undefined southwards   */
            return (row >= ih - 1) ? defaultVal : *(ptr + iw);
        }

        unsigned char GetNorthUChar(unsigned char *ptr, const int row, const int iw, unsigned char defaultVal)
        {
            /* catch case where image is undefined northwards   */
            return (row < 1) ? defaultVal : *(ptr - iw);
        }

        unsigned char GetEastUChar(unsigned char *ptr, const int col, const int iw, unsigned char defaultVal)
        {
            /* catch case where image is undefined eastwards   */
            return (col >= iw - 1) ? defaultVal : *(ptr + 1);
        }

        unsigned char GetWestUChar(unsigned char *ptr, const int col, unsigned char defaultVal)
        {
            /* catch case where image is undefined westwards   */
            return (col < 1) ? defaultVal : *(ptr - 1);
        }

        void RemoveInconDirs(int *imap, const int mw, const int mh,
                    const Dir2Rad& dir2rad, const LQMParams& lqmParams) {

            /* Compute center coords of IMAP */
            int cx = mw>>1;
            int cy = mh>>1;

            /* Do pass, while directions have been removed in a pass ... */
            int nremoved = 0;
            do {
                /* Reinitialize number of removed directions to 0 */
                nremoved = 0;

                /* Start at center */
                int* iptr = imap + (cy * mw) + cx;
                /* If valid IMAP direction and test for removal is true ... */
                if ((*iptr != LQM_INVALID_DIR)&&
                    (RemoveDir(imap, cx, cy, mw, mh, dir2rad, lqmParams))){

                    /* Set to INVALID */
                    *iptr = LQM_INVALID_DIR;
                    /* Bump number of removed IMAP directions */
                    nremoved++;
                }

                /* Initialize side indices of concentric boxes */
                int lbox = cx-1;
                int tbox = cy-1;
                int rbox = cx+1;
                int bbox = cy+1;

                /* Grow concentric boxes, until ALL edges of imap are exceeded */
                while ((lbox >= 0) || (rbox < mw) || (tbox >= 0) || (bbox < mh)) {

                    /* test top edge of box */
                    if (tbox >= 0) {
                        nremoved += TestTopEdge(lbox, tbox, rbox, imap, mw, mh,
                                        dir2rad, lqmParams);
                    }

                    /* test right edge of box */
                    if (rbox < mw) {
                        nremoved += TestRightEdge(tbox, rbox, bbox, imap, mw, mh,
                                        dir2rad, lqmParams);
                    }

                    /* test bottom edge of box */
                    if (bbox < mh) {
                        nremoved += TestBottomEdge(lbox, rbox, bbox, imap, mw, mh,
                                        dir2rad, lqmParams);
                    }

                    /* test left edge of box */
                    if (lbox >=0) {
                        nremoved += TestLeftEdge(lbox, tbox, bbox, imap, mw, mh,
                                        dir2rad, lqmParams);
                    }

                    /* Resize current box */
                    lbox--;
                    tbox--;
                    rbox++;
                    bbox++;
                }
            } while (nremoved);
        }

        int RemoveDir(int *imap, const int mx, const int my,
                    const int mw, const int mh, const Dir2Rad& dir2rad,
                    const LQMParams& lqmParams) {

            int avrdir, nvalid;
            double dir_strength;

            /* Compute average direction from neighbors, returning the */
            /* number of valid neighbors used in the computation, and  */
            /* the "strength" of the average direction.                */
            Average8NeighDir(&avrdir, &dir_strength, &nvalid, imap, mx, my, mw, mh,
                                dir2rad);

            /* Conduct valid neighbor test (Ex. thresh==3) */
            if (nvalid < lqmParams.rmv_valid_nbr_min) {
                return 1;
            }

            /* If strength of average neighbor direction is large enough to */
            /* put credence in ... (Ex. thresh==0.2)                        */
            if (dir_strength >= lqmParams.dir_strength_min) {

                /* Conduct direction distance test (Ex. thresh==3) */
                /* Compute minimum absolute distance between current and       */
                /* average directions accounting for wrapping from 0 to NDIRS. */
                int dist = std::abs(avrdir - *(imap+(my*mw)+mx));
                dist = std::min(dist, dir2rad.ndirs-dist);
                if (dist > lqmParams.dir_distance_max) {
                    return 2;
                }
            }

            /* Otherwise, the strength of the average direciton is not strong enough */
            /* to put credence in, so leave the current block's directon alone.      */

            return 0;
        }

        void Average8NeighDir(int *avrdir, double *dir_strength, int *nvalid,
                            int *imap, const int mx, const int my,
                            const int mw, const int mh,
                            const Dir2Rad& dir2rad) {

            /* Compute neighbor coordinates to current IMAP direction */
            int e = mx+1;  /* East */
            int w = mx-1;  /* West */
            int n = my-1;  /* North */
            int s = my+1;  /* South */

            /* Intialize accumulators */
            *nvalid = 0;
            double cospart = 0.0;
            double sinpart = 0.0;

            /* 1. Test NW */
            /* If NW point within IMAP boudaries ... */
            if ((w >= 0) && (n >= 0)) {
                int* iptr = imap + (n*mw) + w;
                /* If valid direction ... */
                if (*iptr != LQM_INVALID_DIR){
                    /* Accumulate cosine and sine components of the direction */
                    cospart += dir2rad.cosLut[static_cast<std::size_t>(*iptr)];
                    sinpart += dir2rad.sinLut[static_cast<std::size_t>(*iptr)];
                    /* Bump number of accumulated directions */
                    (*nvalid)++;
                }
            }

            /* 2. Test N */
            /* If N point within IMAP boudaries ... */
            if (n >= 0) {
                int* iptr = imap + (n*mw) + mx;
                /* If valid direction ... */
                if (*iptr != LQM_INVALID_DIR){
                    /* Accumulate cosine and sine components of the direction */
                    cospart += dir2rad.cosLut[static_cast<std::size_t>(*iptr)];
                    sinpart += dir2rad.sinLut[static_cast<std::size_t>(*iptr)];
                    /* Bump number of accumulated directions */
                    (*nvalid)++;
                }
            }

            /* 3. Test NE */
            /* If NE point within IMAP boudaries ... */
            if ((e < mw) && (n >= 0)) {
                int* iptr = imap + (n*mw) + e;
                /* If valid direction ... */
                if(*iptr != LQM_INVALID_DIR){
                    /* Accumulate cosine and sine components of the direction */
                    cospart += dir2rad.cosLut[static_cast<std::size_t>(*iptr)];
                    sinpart += dir2rad.sinLut[static_cast<std::size_t>(*iptr)];
                    /* Bump number of accumulated directions */
                    (*nvalid)++;
                }
            }

            /* 4. Test E */
            /* If E point within IMAP boudaries ... */
            if (e < mw) {
                int* iptr = imap + (my*mw) + e;
                /* If valid direction ... */
                if (*iptr != LQM_INVALID_DIR){
                    /* Accumulate cosine and sine components of the direction */
                    cospart += dir2rad.cosLut[static_cast<std::size_t>(*iptr)];
                    sinpart += dir2rad.sinLut[static_cast<std::size_t>(*iptr)];
                    /* Bump number of accumulated directions */
                    (*nvalid)++;
                }
            }

            /* 5. Test SE */
            /* If SE point within IMAP boudaries ... */
            if ((e < mw) && (s < mh)) {
                int* iptr = imap + (s*mw) + e;
                /* If valid direction ... */
                if (*iptr != LQM_INVALID_DIR) {
                    /* Accumulate cosine and sine components of the direction */
                    cospart += dir2rad.cosLut[static_cast<std::size_t>(*iptr)];
                    sinpart += dir2rad.sinLut[static_cast<std::size_t>(*iptr)];
                    /* Bump number of accumulated directions */
                    (*nvalid)++;
                }
            }

            /* 6. Test S */
            /* If S point within IMAP boudaries ... */
            if (s < mh) {
                int* iptr = imap + (s*mw) + mx;
                /* If valid direction ... */
                if (*iptr != LQM_INVALID_DIR) {
                    /* Accumulate cosine and sine components of the direction */
                    cospart += dir2rad.cosLut[static_cast<std::size_t>(*iptr)];
                    sinpart += dir2rad.sinLut[static_cast<std::size_t>(*iptr)];
                    /* Bump number of accumulated directions */
                    (*nvalid)++;
                }
            }

            /* 7. Test SW */
            /* If SW point within IMAP boudaries ... */
            if ((w >= 0) && (s < mh)) {
                int* iptr = imap + (s*mw) + w;
                /* If valid direction ... */
                if (*iptr != LQM_INVALID_DIR) {
                    /* Accumulate cosine and sine components of the direction */
                    cospart += dir2rad.cosLut[static_cast<std::size_t>(*iptr)];
                    sinpart += dir2rad.sinLut[static_cast<std::size_t>(*iptr)];
                    /* Bump number of accumulated directions */
                    (*nvalid)++;
                }
            }

            /* 8. Test W */
            /* If W point within IMAP boudaries ... */
            if (w >= 0) {
                int* iptr = imap + (my*mw) + w;
                /* If valid direction ... */
                if (*iptr != LQM_INVALID_DIR) {
                    /* Accumulate cosine and sine components of the direction */
                    cospart += dir2rad.cosLut[static_cast<std::size_t>(*iptr)];
                    sinpart += dir2rad.sinLut[static_cast<std::size_t>(*iptr)];
                    /* Bump number of accumulated directions */
                    (*nvalid)++;
                }
            }

            /* If there were no neighbors found with valid direction ... */
            if (*nvalid == 0) {
                /* Return INVALID direction. */
                *dir_strength = 0;
                *avrdir = LQM_INVALID_DIR;
                return;
            }

            /* Compute averages of accumulated cosine and sine direction components */
            cospart /= static_cast<double>(*nvalid);
            sinpart /= static_cast<double>(*nvalid);

            /* Compute directional strength as hypotenuse (without sqrt) of average */
            /* cosine and sine direction components.  Believe this value will be on */
            /* the range of [0 .. 1].                                               */
            *dir_strength = (cospart * cospart) + (sinpart * sinpart);
            /* Need to truncate precision so that answers are consistent   */
            /* on different computer architectures when comparing doubles. */
            *dir_strength = trunc_dbl_precision(*dir_strength);

            /* If the direction strength is not sufficiently high ... */
            if (*dir_strength < LQM_DIR_STRENGTH_MIN) {
                /* Return INVALID direction. */
                *dir_strength = 0;
                *avrdir = LQM_INVALID_DIR;
                return;
            }

            /* Compute angle (in radians) from Arctan of avarage         */
            /* cosine and sine direction components.  I think this order */
            /* is necessary because 0 direction is vertical and positive */
            /* direction is clockwise.                                   */
            double theta = atan2(sinpart, cospart);

            /* Atan2 returns theta on range [-PI..PI].  Adjust theta so that */
            /* it is on the range [0..2PI].                                  */
            double pi2 = 2*LQM_PI;
            theta += pi2;
            theta = fmod(theta, pi2);

            /* Pi_factor sets the period of the trig functions to NDIRS units in x. */
            /* For example, if NDIRS==16, then pi_factor = 2(PI/16) = .3926...      */
            /* Dividing theta (in radians) by this factor ((1/pi_factor)==2.546...) */
            /* will produce directions on the range [0..NDIRS].                     */
            double pi_factor = pi2/static_cast<double>(dir2rad.ndirs); /* 2(M_PI/ndirs) */

            /* Round off the direction and return it as an average direction */
            /* for the neighborhood.                                         */
            double avr = theta / pi_factor;
            /* Need to truncate precision so that answers are consistent */
            /* on different computer architectures when rounding doubles. */
            avr = trunc_dbl_precision(avr);
            *avrdir = sround(avr);

            /* Really do need to map values > NDIRS back onto [0..NDIRS) range. */
            *avrdir %= dir2rad.ndirs;
        }

        int TestTopEdge(const int lbox, const int tbox, const int rbox,
                        int *imap, const int mw, const int mh,
                        const Dir2Rad& dir2rad, const LQMParams& lqmParams) {

            /* Initialize number of directions removed on edge to 0 */
            int nremoved = 0;

            /* Set start pointer to top-leftmost point of box, or set it to */
            /* the leftmost point in the IMAP row (0), whichever is larger. */
            int sx = std::max(lbox, 0);
            int* sptr = imap + (tbox*mw) + sx;

            /* Set end pointer to either 1 point short of the top-rightmost */
            /* point of box, or set it to the rightmost point in the IMAP   */
            /* row (lastx=mw-1), whichever is smaller.                      */
            int ex = std::min(rbox-1, mw-1);
            int* eptr = imap + (tbox*mw) + ex;

            /* For each point on box's edge ... */
            for (int *iptr = sptr, bx = sx, by = tbox;
                iptr <= eptr;
                iptr++, bx++) {
                /* If valid IMAP direction and test for removal is true ... */
                if ((*iptr != LQM_INVALID_DIR) &&
                    (RemoveDir(imap, bx, by, mw, mh, dir2rad, lqmParams))) {
                    /* Set to INVALID */
                    *iptr = LQM_INVALID_DIR;
                    /* Bump number of removed IMAP directions */
                    nremoved++;
                }
            }

            /* Return the number of directions removed on edge */
            return nremoved;
        }

        int TestRightEdge(const int tbox, const int rbox,
                            const int bbox, int *imap, const int mw, const int mh,
                            const Dir2Rad& dir2rad, const LQMParams& lqmParams) {

            /* Initialize number of directions removed on edge to 0 */
            int nremoved = 0;

            /* Set start pointer to top-rightmost point of box, or set it to */
            /* the topmost point in IMAP column (0), whichever is larger.    */
            int sy = std::max(tbox, 0);
            int* sptr = imap + (sy*mw) + rbox;

            /* Set end pointer to either 1 point short of the bottom-    */
            /* rightmost point of box, or set it to the bottommost point */
            /* in the IMAP column (lasty=mh-1), whichever is smaller.    */
            int ey = std::min(bbox-1,mh-1);
            int* eptr = imap + (ey*mw) + rbox;

            /* For each point on box's edge ... */
            for (int *iptr = sptr, bx = rbox, by = sy;
                iptr <= eptr;
                iptr+=mw, by++) {
                /* If valid IMAP direction and test for removal is true ... */
                if ((*iptr != LQM_INVALID_DIR)&&
                    (RemoveDir(imap, bx, by, mw, mh, dir2rad, lqmParams))) {
                    /* Set to INVALID */
                    *iptr = LQM_INVALID_DIR;
                    /* Bump number of removed IMAP directions */
                    nremoved++;
                }
            }

            /* Return the number of directions removed on edge */
            return nremoved;
        }

        int TestBottomEdge(const int lbox, const int rbox,
                            const int bbox, int *imap, const int mw, const int mh,
                            const Dir2Rad& dir2rad, const LQMParams& lqmParams) {

            /* Initialize number of directions removed on edge to 0 */
            int nremoved = 0;

            /* Set start pointer to bottom-rightmost point of box, or set it to the */
            /* rightmost point in the IMAP ROW (lastx=mw-1), whichever is smaller.  */
            int sx = std::min(rbox, mw-1);
            int* sptr = imap + (bbox*mw) + sx;

            /* Set end pointer to either 1 point short of the bottom-    */
            /* lefttmost point of box, or set it to the leftmost point   */
            /* in the IMAP row (x=0), whichever is larger.               */
            int ex = std::max(lbox-1, 0);
            int* eptr = imap + (bbox*mw) + ex;

            /* For each point on box's edge ... */
            for (int *iptr = sptr, bx = sx, by = bbox;
                iptr >= eptr;
                iptr--, bx--) {
                /* If valid IMAP direction and test for removal is true ... */
                if ((*iptr != LQM_INVALID_DIR)&&
                    (RemoveDir(imap, bx, by, mw, mh, dir2rad, lqmParams))) {
                    /* Set to INVALID */
                    *iptr = LQM_INVALID_DIR;
                    /* Bump number of removed IMAP directions */
                    nremoved++;
                }
            }

            /* Return the number of directions removed on edge */
            return nremoved;
        }

        int TestLeftEdge(const int lbox, const int tbox,
                        const int bbox, int *imap, const int mw, const int mh,
                        const Dir2Rad& dir2rad, const LQMParams& lqmParams) {

            /* Initialize number of directions removed on edge to 0 */
            int nremoved = 0;

            /* Set start pointer to bottom-leftmost point of box, or set it to */
            /* the bottommost point in IMAP column (lasty=mh-1), whichever     */
            /* is smaller.                                                     */
            int sy = std::min(bbox, mh-1);
            int* sptr = imap + (sy*mw) + lbox;

            /* Set end pointer to either 1 point short of the top-leftmost */
            /* point of box, or set it to the topmost point in the IMAP    */
            /* column (y=0), whichever is larger.                          */
            int ey = std::max(tbox-1, 0);
            int* eptr = imap + (ey*mw) + lbox;

            /* For each point on box's edge ... */
            for (int *iptr = sptr, bx = lbox, by = sy;
                iptr >= eptr;
                iptr-=mw, by--) {
                /* If valid IMAP direction and test for removal is true ... */
                if ((*iptr != LQM_INVALID_DIR)&&
                    (RemoveDir(imap, bx, by, mw, mh, dir2rad, lqmParams))) {
                    /* Set to INVALID */
                    *iptr = LQM_INVALID_DIR;
                    /* Bump number of removed IMAP directions */
                    nremoved++;
                }
            }

            /* Return the number of directions removed on edge */
            return nremoved;
        }

        void SmoothDirectionMap(int *direction_map, int *low_contrast_map,
                        const int mw, const int mh,
                        const Dir2Rad& dir2rad, const LQMParams& lqmParams) {

            /* Assign pointers to beginning of both maps. */
            int* dptr = direction_map;
            int* cptr = low_contrast_map;

            /* Foreach block in maps ... */
            for (int my = 0; my < mh; ++my) {
                for (int mx = 0; mx < mw; ++mx) {
                    /* If the current block does NOT have LOW CONTRAST ... */
                    if (!*cptr) {

                        /* Compute average direction from neighbors, returning the */
                        /* number of valid neighbors used in the computation, and  */
                        /* the "strength" of the average direction.                */
                        int avrdir;
                        double dir_strength;
                        int nvalid;
                        Average8NeighDir(&avrdir, &dir_strength, &nvalid,
                                        direction_map, mx, my, mw, mh, dir2rad);

                        /* If average direction strength is strong enough */
                        /*    (Ex. thresh==0.2)...                        */
                        if(dir_strength >= lqmParams.dir_strength_min) {
                            /* If Direction Map direction is valid ... */
                            if (*dptr != LQM_INVALID_DIR) {
                                /* Conduct valid neighbor test (Ex. thresh==3)... */
                                if(nvalid >= lqmParams.rmv_valid_nbr_min){
                                    /* Reassign valid direction with average direction. */
                                    *dptr = avrdir;
                                }
                            }
                            /* Otherwise direction is invalid ... */
                            else {
                                /* Even if DIRECTION_MAP value is invalid, if number of */
                                /* valid neighbors is big enough (Ex. thresh==7)...     */
                                if (nvalid >= lqmParams.smth_valid_nbr_min) {
                                    /* Assign invalid direction with average direction. */
                                    *dptr = avrdir;
                                }
                            }
                        }
                    }
                    /* Otherwise, block has LOW CONTRAST, so keep INVALID direction. */

                    /* Bump to next block in maps. */
                    dptr++;
                    cptr++;
                }
            }
        }

        void InterpolateDirectionMap(int *direction_map, int *low_contrast_map,
                            const int mw, const int mh, const LQMParams& lqmParams)
        {
            /* Allocate output (interpolated) Direction Map. */
            std::vector<int> interpolated_dmap(static_cast<std::size_t>(mw * mh));
            int* omap = interpolated_dmap.data();

            int* dptr = direction_map;
            int* cptr = low_contrast_map;
            int* optr = omap;

            /* Foreach block in the maps ... */
            for (int y = 0; y < mh; ++y) {
                for (int x = 0; x < mw; ++x) {
                    /* If image block is NOT LOW CONTRAST and has LQM_INVALID_DIR direction ... */
                    if ((!*cptr) && (*dptr == LQM_INVALID_DIR)) {
                        /* Set neighbor accumulators to 0. */
                        int total_found = 0;
                        int total_dist = 0;

                        int nbr_x;
                        int nbr_y;

                        /* Find north neighbor. */
                        int n_dir; // Only set when n_found, and only accessed later when n_found
                        int n_dist = 0;
                        int n_found = FindValidBlock(&n_dir, &nbr_x, &nbr_y,
                                                        direction_map, low_contrast_map,
                                                        x, y, mw, mh, 0, -1);
                        if (n_found == LQM_FOUND) {
                            /* Compute north distance. */
                            n_dist = y - nbr_y;
                            /* Accumulate neighbor distance. */
                            total_dist += n_dist;
                            /* Bump number of neighbors found. */
                            total_found++;
                        }

                        /* Find east neighbor. */
                        int e_dir;
                        int e_dist = 0;
                        int e_found = FindValidBlock(&e_dir, &nbr_x, &nbr_y,
                                                        direction_map, low_contrast_map,
                                                        x, y, mw, mh, 1, 0);
                        if (e_found == LQM_FOUND) {
                            /* Compute east distance. */
                            e_dist = nbr_x - x;
                            /* Accumulate neighbor distance. */
                            total_dist += e_dist;
                            /* Bump number of neighbors found. */
                            total_found++;
                        }

                        /* Find south neighbor. */
                        int s_dir;
                        int s_dist = 0;
                        int s_found = FindValidBlock(&s_dir, &nbr_x, &nbr_y,
                                                        direction_map, low_contrast_map,
                                                        x, y, mw, mh, 0, 1);
                        if (s_found == LQM_FOUND) {
                            /* Compute south distance. */
                            s_dist = nbr_y - y;
                            /* Accumulate neighbor distance. */
                            total_dist += s_dist;
                            /* Bump number of neighbors found. */
                            total_found++;
                        }

                        /* Find west neighbor. */
                        int w_dir;
                        int w_dist = 0;
                        int w_found = FindValidBlock(&w_dir, &nbr_x, &nbr_y,
                                                        direction_map, low_contrast_map,
                                                        x, y, mw, mh, -1, 0);
                        if (w_found == LQM_FOUND) {
                            /* Compute west distance. */
                            w_dist = x - nbr_x;
                            /* Accumulate neighbor distance. */
                            total_dist += w_dist;
                            /* Bump number of neighbors found. */
                            total_found++;
                        }

                        /* If a sufficient number of neighbors found (Ex. 2) ... */
                        if (total_found >= lqmParams.min_interpolate_nbrs) {
                            /* Accumulate weighted sum of neighboring directions     */
                            /* inversely related to the distance from current block. */
                            int total_delta = 0;

                            int n_delta = 0;
                            int e_delta = 0;
                            int s_delta = 0;
                            int w_delta = 0;
                            /* If neighbor found to the north ... */
                            if (n_found) {
                                n_delta = total_dist - n_dist;
                                total_delta += n_delta;
                            }
                            /* If neighbor found to the east ... */
                            if (e_found) {
                                e_delta = total_dist - e_dist;
                                total_delta += e_delta;
                            }
                            /* If neighbor found to the south ... */
                            if (s_found) {
                                s_delta = total_dist - s_dist;
                                total_delta += s_delta;
                            }
                            /* If neighbor found to the west ... */
                            if (w_found) {
                                w_delta = total_dist - w_dist;
                                total_delta += w_delta;
                            }

                            double avr_dir = 0.0;

                            if(n_found){
                                avr_dir += (n_dir*(n_delta/static_cast<double>(total_delta)));
                            }
                            if(e_found){
                                avr_dir += (e_dir*(e_delta/static_cast<double>(total_delta)));
                            }
                            if(s_found){
                                avr_dir += (s_dir*(s_delta/static_cast<double>(total_delta)));
                            }
                            if(w_found){
                                avr_dir += (w_dir*(w_delta/static_cast<double>(total_delta)));
                            }

                            /* Need to truncate precision so that answers are consistent  */
                            /* on different computer architectures when rounding doubles. */
                            avr_dir = trunc_dbl_precision(avr_dir);

                            // Flag areas with interpolated directions
                            /* Assign interpolated direction to output Direction Map. */
                            int new_dir = sround(avr_dir);

                            *optr = new_dir;
                        }
                        else{
                        /* Otherwise, the direction remains LQM_INVALID. */
                        *optr = *dptr;
                        }
                    }
                    else{
                        /* Otherwise, assign the current direction to the output block. */
                        *optr = *dptr;
                    }

                    /* Bump to the next block in the maps ... */
                    dptr++;
                    cptr++;
                    optr++;
                }
            }

            /* Copy the interpolated directions into the input map. */
            memcpy(direction_map, omap, static_cast<std::size_t>(mw * mh) * sizeof(int));
        }

        int FindValidBlock(int *nbr_dir, int *nbr_x, int *nbr_y,
                            int *direction_map, int *low_contrast_map,
                            const int sx, const int sy,
                            const int mw, const int mh,
                            const int x_incr, const int y_incr)
        {
            /* Initialize starting block coords. */
            int x = sx + x_incr;
            int y = sy + y_incr;

            /* While we are not outside the boundaries of the map ... */
            while ((x >= 0) && (x < mw) && (y >= 0) && (y < mh)) {
                /* Stop unsuccessfully if we encounter a LOW CONTRAST block. */
                // Could change from boolean to continuous with const threshold
                if (*(low_contrast_map+(y*mw)+x)) {
                    return LQM_NOT_FOUND;
                }

                /* Stop successfully if we encounter a block with valid direction. */
                int dir = *(direction_map+(y*mw)+x);
                if (dir >= 0) {
                    *nbr_dir = dir;
                    *nbr_x = x;
                    *nbr_y = y;
                    return LQM_FOUND;
                }

                /* Otherwise, advance to the next block in the map. */
                x += x_incr;
                y += y_incr;
            }

            /* If we get here, then we did not find a valid block in the given */
            /* direction in the map.                                           */
            return LQM_NOT_FOUND;
        }

        void SetMarginBlocks(int *map, const int mw, const int mh,
                            const int margin_value)
        {
            int* ptr1 = map;
            int* ptr2 = map+((mh-1)*mw);
            for (int x = 0; x < mw; ++x) {
                *ptr1++ = margin_value;
                *ptr2++ = margin_value;
            }

            ptr1 = map + mw;
            ptr2 = map + mw + mw - 1;
            for (int y = 1; y < mh-1; ++y) {
                *ptr1 = margin_value;
                *ptr2 = margin_value;
                ptr1 += mw;
                ptr2 += mw;
            }
        }

        void GenHighCurveMap(std::vector<int>& out_hcmap, int *direction_map,
                        const int mw, const int mh, const LQMParams& lqmParams,
                        unsigned char *validNeighbors, unsigned char *curvatureMap)
        {
            int mapsize = mw*mh;

            /* Allocate High Curvature Map initialized to LQM_FALSE (0) */
            out_hcmap = std::vector<int>(static_cast<std::size_t>(mapsize), LQM_FALSE);

            int* hptr = out_hcmap.data();
            int* dptr = direction_map;

            /* Foreach row in maps ... */
            for (int by = 0; by < mh; ++by) {
                /* Foreach column in maps ... */
                for(int bx = 0; bx < mw; ++bx) {
                    int offset = by*mw + bx;

                    curvatureMap[offset]=255;

                    /* Count number of valid neighbors around current block ... */
                    int nvalid = NumValid8Neigh(direction_map, bx, by, mw, mh);
                    validNeighbors[offset]= static_cast<unsigned char>(nvalid * 31); //rescale 0..8 to 0..248

                    /* If valid neighbors exist ... */
                    if (nvalid > 0) {
                        /* If current block's direction is INVALID ... */
                        if (*dptr == LQM_INVALID_DIR) {
                            /* If a sufficient number of VALID neighbors exists ... */
                            if (nvalid >= lqmParams.vort_valid_nbr_min) {
                                /* Measure vorticity of neighbors. */
                                int vmeasure = vorticity(direction_map, bx, by, mw, mh,
                                                    lqmParams.num_directions);

                                /* If vorticity is sufficiently high ... */
                                if (vmeasure >= lqmParams.highcurv_vorticity_min) {
                                    /* Flag block as HIGH CURVATURE. */
                                    *hptr = LQM_TRUE;
                                }
                            }
                        }
                        /* Otherwise block has valid direction ... */
                        else {
                            /* Measure curvature around the valid block. */
                            int cmeasure = curvature(direction_map, bx, by, mw, mh,
                                                    lqmParams.num_directions);

                            curvatureMap[static_cast<std::size_t>(offset)]= static_cast<unsigned char> (std::min(254,std::max(0,cmeasure*16)));  //255=undefined; *16 to rescale

                            /* If curvature is sufficiently high ... */
                            if (cmeasure >= lqmParams.highcurv_curvature_min) {
                                *hptr = LQM_TRUE;
                            }
                        }
                    } /* Else (nvalid <= 0) */

                    /* Bump pointers to next block in maps. */
                    dptr++;
                    hptr++;

                } /* bx */
            } /* by */
        }

        int NumValid8Neigh(int *imap, const int mx, const int my,
                            const int mw, const int mh)
        {
            /* Initialize VALID IMAP counter to zero. */
            int nvalid = 0;

            /* Compute neighbor coordinates to current IMAP direction */
            int e_ind = mx+1;  /* East index */
            int w_ind = mx-1;  /* West index */
            int n_ind = my-1;  /* North index */
            int s_ind = my+1;  /* South index */

            /* 1. Test NW IMAP value.  */
            /* If neighbor indices are within IMAP boundaries and it is VALID ... */
            if ((w_ind >= 0) && (n_ind >= 0) && (*(imap + (n_ind*mw) + w_ind) >= 0)) {
                /* Bump VALID counter. */
                nvalid++;
            }

            /* 2. Test N IMAP value.  */
            if((n_ind >= 0) && (*(imap + (n_ind*mw) + mx) >= 0)) {
                nvalid++;
            }

            /* 3. Test NE IMAP value. */
            if ((n_ind >= 0) && (e_ind < mw) && (*(imap + (n_ind*mw) + e_ind) >= 0)) {
                nvalid++;
            }

            /* 4. Test E IMAP value. */
            if ((e_ind < mw) && (*(imap + (my*mw) + e_ind) >= 0)) {
                nvalid++;
            }

            /* 5. Test SE IMAP value. */
            if ((e_ind < mw) && (s_ind < mh) && (*(imap + (s_ind*mw) + e_ind) >= 0)) {
                nvalid++;
            }

            /* 6. Test S IMAP value. */
            if ((s_ind < mh) && (*(imap + (s_ind*mw) + mx) >= 0)) {
                nvalid++;
            }

            /* 7. Test SW IMAP value. */
            if((w_ind >= 0) && (s_ind < mh) && (*(imap + (s_ind*mw) + w_ind) >= 0))
                nvalid++;

            /* 8. Test W IMAP value. */
            if ((w_ind >= 0) && (*(imap + (my*mw) + w_ind) >= 0)) {
                nvalid++;
            }

            /* Return number of neighbors with VALID IMAP values. */
            return nvalid;
        }

        int vorticity(int *imap, const int mx, const int my,
                    const int mw, const int mh, const int ndirs)
        {
            /* Compute neighbor coordinates to current IMAP direction */
            int e_ind = mx+1;  /* East index */
            int w_ind = mx-1;  /* West index */
            int n_ind = my-1;  /* North index */
            int s_ind = my+1;  /* South index */

            int nw_val;
            int n_val;
            int ne_val;
            int e_val;
            int se_val;
            int s_val;
            int sw_val;
            int w_val;

            /* 1. Get NW IMAP value.  */
            /* If neighbor indices are within IMAP boundaries ... */
            if ((w_ind >= 0) && (n_ind >= 0)) {
                /* Set neighbor value to IMAP value. */
                nw_val = *(imap + (n_ind*mw) + w_ind);
            }
            else {
                /* Otherwise, set the neighbor value to INVALID. */
                nw_val = LQM_INVALID_DIR;
            }

            /* 2. Get N IMAP value.  */
            if (n_ind >= 0) {
                n_val = *(imap + (n_ind*mw) + mx);
            }
            else {
                n_val = LQM_INVALID_DIR;
            }

            /* 3. Get NE IMAP value. */
            if ((n_ind >= 0) && (e_ind < mw)) {
                ne_val = *(imap + (n_ind*mw) + e_ind);
            }
            else {
                ne_val = LQM_INVALID_DIR;
            }

            /* 4. Get E IMAP value. */
            if (e_ind < mw) {
                e_val = *(imap + (my*mw) + e_ind);
            }
            else {
                e_val = LQM_INVALID_DIR;
            }

            /* 5. Get SE IMAP value. */
            if ((e_ind < mw) && (s_ind < mh)) {
                se_val = *(imap + (s_ind*mw) + e_ind);
            }
            else {
                se_val = LQM_INVALID_DIR;
            }

            /* 6. Get S IMAP value. */
            if(s_ind < mh) {
                s_val = *(imap + (s_ind*mw) + mx);
            }
            else {
                s_val = LQM_INVALID_DIR;
            }

            /* 7. Get SW IMAP value. */
            if ((w_ind >= 0) && (s_ind < mh)) {
                sw_val = *(imap + (s_ind*mw) + w_ind);
            }
            else {
                sw_val = LQM_INVALID_DIR;
            }

            /* 8. Get W IMAP value. */
            if (w_ind >= 0) {
                w_val = *(imap + (my*mw) + w_ind);
            }
            else {
                w_val = LQM_INVALID_DIR;
            }

            /* Now that we have all IMAP neighbors, accumulate vorticity between */
            /* the neighboring directions.                                       */

            /* Initialize vorticity accumulator to zero. */
            int vmeasure = 0;

            /* 1. NW & N */
            AccumNeighVorticity(&vmeasure, nw_val, n_val, ndirs);

            /* 2. N & NE */
            AccumNeighVorticity(&vmeasure, n_val, ne_val, ndirs);

            /* 3. NE & E */
            AccumNeighVorticity(&vmeasure, ne_val, e_val, ndirs);

            /* 4. E & SE */
            AccumNeighVorticity(&vmeasure, e_val, se_val, ndirs);

            /* 5. SE & S */
            AccumNeighVorticity(&vmeasure, se_val, s_val, ndirs);

            /* 6. S & SW */
            AccumNeighVorticity(&vmeasure, s_val, sw_val, ndirs);

            /* 7. SW & W */
            AccumNeighVorticity(&vmeasure, sw_val, w_val, ndirs);

            /* 8. W & NW */
            AccumNeighVorticity(&vmeasure, w_val, nw_val, ndirs);

            /* Return the accumulated vorticity measure. */
            return vmeasure;
        }

        void AccumNeighVorticity(int *vmeasure, const int dir1, const int dir2,
                                const int ndirs)
        {
            /* Measure difference in direction between a pair of neighboring */
            /* directions.                                                   */
            /* If both neighbors are not equal and both are VALID ... */
            if ((dir1 != dir2) && (dir1 >= 0)&&(dir2 >= 0)) {
                /* Measure the clockwise distance from the first to the second */
                /* directions.                                                 */
                int dist = dir2 - dir1;
                /* If dist is negative, then clockwise distance must wrap around */
                /* the high end of the direction range. For example:             */
                /*              dir1 = 8                                         */
                /*              dir2 = 3                                         */
                /*       and   ndirs = 16                                        */
                /*             3 - 8 = -5                                        */
                /*        so  16 - 5 = 11  (the clockwise distance from 8 to 3)  */
                if (dist < 0) {
                    dist += ndirs;
                }
                /* If the change in clockwise direction is larger than 90 degrees as */
                /* in total the total number of directions covers 180 degrees.       */
                if (dist > (ndirs>>1)) {
                    /* Decrement the vorticity measure. */
                    (*vmeasure)--;
                }
                else {
                    /* Otherwise, bump the vorticity measure. */
                    (*vmeasure)++;
                }
            }
            /* Otherwise both directions are either equal or  */
            /* one or both directions are INVALID, so ignore. */
        }

        int curvature(int *imap, const int mx, const int my,
                    const int mw, const int mh, const int ndirs)
        {
            int nw_val, n_val, ne_val, e_val, se_val, s_val, sw_val, w_val;

            /* Compute neighbor coordinates to current IMAP direction */
            int e_ind = mx+1;  /* East index */
            int w_ind = mx-1;  /* West index */
            int n_ind = my-1;  /* North index */
            int s_ind = my+1;  /* South index */

            /* 1. Get NW IMAP value.  */
            /* If neighbor indices are within IMAP boundaries ... */
            if ((w_ind >= 0) && (n_ind >= 0)) {
                /* Set neighbor value to IMAP value. */
                nw_val = *(imap + (n_ind*mw) + w_ind);
            }
            else {
                /* Otherwise, set the neighbor value to INVALID. */
                nw_val = LQM_INVALID_DIR;
            }

            /* 2. Get N IMAP value.  */
            if (n_ind >= 0) {
                n_val = *(imap + (n_ind*mw) + mx);
            }
            else {
                n_val = LQM_INVALID_DIR;
            }

            /* 3. Get NE IMAP value. */
            if ((n_ind >= 0) && (e_ind < mw)) {
                ne_val = *(imap + (n_ind*mw) + e_ind);
            }
            else {
                ne_val = LQM_INVALID_DIR;
            }

            /* 4. Get E IMAP value. */
            if (e_ind < mw) {
                e_val = *(imap + (my*mw) + e_ind);
            }
            else {
                e_val = LQM_INVALID_DIR;
            }

            /* 5. Get SE IMAP value. */
            if ((e_ind < mw) && (s_ind < mh)) {
                se_val = *(imap + (s_ind*mw) + e_ind);
            }
            else {
                se_val = LQM_INVALID_DIR;
            }

            /* 6. Get S IMAP value. */
            if (s_ind < mh) {
                s_val = *(imap + (s_ind*mw) + mx);
            }
            else {
                s_val = LQM_INVALID_DIR;
            }

            /* 7. Get SW IMAP value. */
            if ((w_ind >= 0) && (s_ind < mh)) {
                sw_val = *(imap + (s_ind*mw) + w_ind);
            }
            else {
                sw_val = LQM_INVALID_DIR;
            }

            /* 8. Get W IMAP value. */
            if (w_ind >= 0) {
                w_val = *(imap + (my*mw) + w_ind);
            }
            else {
                w_val = LQM_INVALID_DIR;
            }

            /* Now that we have all IMAP neighbors, determine largest change in */
            /* direction from current block to each of its 8 VALID neighbors.   */

            /* Initialize pointer to current IMAP value. */
            int* iptr = imap + (my*mw) + mx;

            /* Initialize curvature measure to negative as ClosestDirDist() */
            /* always returns -1=LQM_INVALID_DIR or a positive value.                 */
            int cmeasure = -1;

            /* 1. With NW */
            /* Compute closest distance between neighboring directions. */
            int dist = ClosestDirDist(*iptr, nw_val, ndirs);
            /* Keep track of maximum. */
            if (dist > cmeasure) {
                cmeasure = dist;
            }

            /* 2. With N */
            dist = ClosestDirDist(*iptr, n_val, ndirs);
            if (dist > cmeasure) {
                cmeasure = dist;
            }

            /* 3. With NE */
            dist = ClosestDirDist(*iptr, ne_val, ndirs);
            if (dist > cmeasure) {
                cmeasure = dist;
            }

            /* 4. With E */
            dist = ClosestDirDist(*iptr, e_val, ndirs);
            if (dist > cmeasure) {
                cmeasure = dist;
            }

            /* 5. With SE */
            dist = ClosestDirDist(*iptr, se_val, ndirs);
            if (dist > cmeasure) {
                cmeasure = dist;
            }

            /* 6. With S */
            dist = ClosestDirDist(*iptr, s_val, ndirs);
            if (dist > cmeasure) {
                cmeasure = dist;
            }

            /* 7. With SW */
            dist = ClosestDirDist(*iptr, sw_val, ndirs);
            if (dist > cmeasure) {
                cmeasure = dist;
            }

            /* 8. With W */
            dist = ClosestDirDist(*iptr, w_val, ndirs);
            if (dist > cmeasure) {
                cmeasure = dist;
            }

            /* Return maximum difference between current block's IMAP direction */
            /* and the rest of its VALID neighbors.                             */
            return cmeasure;
        }

        int ClosestDirDist(const int dir1, const int dir2, const int ndirs)
        {
            /* Initialize distance to -1 = INVALID. */
            int dist = LQM_INVALID_DIR;

            /* Measure shortest distance between to directions. */
            /* If both neighbors are VALID ... */
            if ((dir1 >= 0) && (dir2 >= 0)) {
                /* Compute inner and outer distances to account for distances */
                /* that wrap around the end of the range of directions, which */
                /* may in fact be closer.                                     */
                int d1 = std::abs(dir2 - dir1);
                int d2 = ndirs - d1;
                dist = std::min(d1, d2);
            }
            /* Otherwise one or both directions are INVALID, so ignore */
            /* and return LQM_INVALID_DIR. */

            /* Return determined closest distance. */
            return dist;
        }

        void GenDirectionChangeMap(unsigned char *directionChangeMap, int *initialDirectionMap, int *direction_map, int mw, int mh)
        {
            for (int counter=0; counter < mw * mh; ++counter) {
                if (initialDirectionMap[counter] < 0 || direction_map[counter] < 0) {
                    directionChangeMap[counter] = 255;
                }
                else {
                    int tmp = std::abs(initialDirectionMap[counter] - direction_map[counter]);
                    directionChangeMap[counter] = static_cast<unsigned char>(16 * std::min(15 - tmp, tmp)); //absolute change in 16 directions: 0-7 - rescaled for visibility in byte
                }
            }
        }

        void Binarize(std::vector<unsigned char>& odata, int& ow, int& oh,
                unsigned char *pdata, const int pw, const int ph,
                int *direction_map, const int mw,
                const RotGrids& dirbingrids, const LQMParams& lqmParams)
        {
            /* 1. Binarize the padded input image using directional block info. */
            BinarizeImage(
                odata, ow, oh, pdata, pw, ph,
                direction_map, mw,
                lqmParams.blocksize, dirbingrids
            );

            /* 2. Fill black and white holes in binary image. */
            /* Scan the binary image, filling holes, 3 times. */
            for (int i = 0; i < lqmParams.num_fill_holes; ++i) {
                FillHoles(odata.data(), ow, oh);
            }
        }

        void BinarizeImage(std::vector<unsigned char>& odata, int& ow, int& oh,
                        unsigned char *pdata, const int pw, const int ph,
                        const int *direction_map, const int mw,
                        const int blocksize, const RotGrids& dirbingrids)
        {
            /* Compute dimensions of "unpadded" binary image results. */
            ow = pw - (dirbingrids.pad<<1);
            oh = ph - (dirbingrids.pad<<1);

            odata.reserve(static_cast<std::size_t>(ow * oh));

            unsigned char* spptr = pdata + (dirbingrids.pad * pw) + dirbingrids.pad;
            for (int iy = 0; iy < oh; iy++) {
                /* Set pixel pointer to start of next row in grid. */
                unsigned char* pptr = spptr;
                for(int ix = 0; ix < ow; ++ix){

                    /* Compute which block the current pixel is in. */
                    int bx = static_cast<int>(ix/blocksize);
                    int by = static_cast<int>(iy/blocksize);
                    /* Get corresponding value in Direction Map. */
                    int mapval = *(direction_map + (by*mw) + bx);

                    unsigned char bVal;
                    /* If current block has has INVALID direction ... */
                    if (mapval == LQM_INVALID_DIR) {
                        /* Set binary pixel to white (255). */
                        bVal = LQM_WHITE_PIXEL;
                    }
                    /* Otherwise, if block has a valid direction ... */
                    else {
                        /* Use directional binarization based on block's direction. */
                        bVal = static_cast<unsigned char>(BinarizeFromDirection(pptr, mapval, dirbingrids));
                    }

                    odata.push_back(bVal);

                    /* Bump input pixel pointer. */
                    pptr++;
                }
                /* Bump pointer to the next row in padded input image. */
                spptr += pw;
            }
        }

        int BinarizeFromDirection(const unsigned char *pptr, const int idir,
                        const RotGrids& dirbingrids)
        {
            int csum = 0;

            /* Assign nickname pointer. */
            const int* grid = dirbingrids.grids[static_cast<std::size_t>(idir)].data();
            /* Calculate center (0-oriented) row in grid. */
            double dcy = (dirbingrids.grid_h-1)/static_cast<double>(2.0);
            /* Need to truncate precision so that answers are consistent */
            /* on different computer architectures when rounding doubles. */
            dcy = trunc_dbl_precision(dcy);
            int cy = sround(dcy);
            /* Initialize grid's pixel offset index to zero. */
            int gi = 0;
            /* Initialize grid's pixel accumulator to zero */
            int gsum = 0;

            /* Foreach row in grid ... */
            for (int gy = 0; gy < dirbingrids.grid_h; ++gy) {
                /* Initialize row pixel sum to zero. */
                int rsum = 0;
                /* Foreach column in grid ... */
                for (int gx = 0; gx < dirbingrids.grid_w; ++gx) {
                    /* Accumulate next pixel along rotated row in grid. */
                    rsum += *(pptr+grid[gi]);
                    /* Bump grid's pixel offset index. */
                    gi++;
                }
                /* Accumulate row sum into grid pixel sum. */
                gsum += rsum;
                /* If current row is center row, then save row sum separately. */
                if (gy == cy) {
                    csum = rsum;
                }
            }

            /* If the center row sum treated as an average is less than the */
            /* total pixel sum in the rotated grid ...                      */
            if ((csum * dirbingrids.grid_h) < gsum) {
                /* Set the binary pixel to BLACK. */
                return LQM_BLACK_PIXEL;
            }
            else {
                /* Otherwise set the binary pixel to WHITE. */
                return LQM_WHITE_PIXEL;
            }
        }

        void FillHoles(unsigned char *bdata, const int iw, const int ih)
        {
            /* 1. Fill 1-pixel wide holes in horizontal runs first ... */
            unsigned char* sptr = bdata + 1;
            /* Foreach row in image ... */
            for (int iy = 0; iy < ih; iy++) {
                /* Initialize pointers to start of next line ... */
                unsigned char* lptr = sptr-1;   /* Left pixel   */
                unsigned char* mptr = sptr;     /* Middle pixel */
                unsigned char* rptr = sptr+1;   /* Right pixel  */
                /* Foreach column in image (less far left and right pixels) ... */
                for(int ix = 1; ix < iw-1; ix++){
                    /* Do we have a horizontal hole of length 1? */
                    if((*lptr != *mptr) && (*lptr == *rptr)){
                        /* If so, then fill it. */
                        *mptr = *lptr;
                        /* Bump passed right pixel because we know it will not */
                        /* be a hole.                                          */
                        lptr+=2;
                        mptr+=2;
                        rptr+=2;
                        /* We bump ix once here and then the FOR bumps it again. */
                        ix++;
                    }
                    else{
                        /* Otherwise, bump to the next pixel to the right. */
                        lptr++;
                        mptr++;
                        rptr++;
                    }
                }
                /* Bump to start of next row. */
                sptr += iw;
            }

            /* 2. Now, fill 1-pixel wide holes in vertical runs ... */
            int iw2 = iw<<1;
            /* Start processing column one row down from the top of the image. */
            sptr = bdata + iw;
            /* Foreach column in image ... */
            for(int ix = 0; ix < iw; ix++){
                /* Initialize pointers to start of next column ... */
                unsigned char* tptr = sptr-iw;   /* Top pixel     */
                unsigned char* mptr = sptr;      /* Middle pixel  */
                unsigned char* bptr = sptr+iw;   /* Bottom pixel  */
                /* Foreach row in image (less top and bottom row) ... */
                for(int iy = 1; iy < ih-1; iy++){
                    /* Do we have a vertical hole of length 1? */
                    if((*tptr != *mptr) && (*tptr == *bptr)){
                        /* If so, then fill it. */
                        *mptr = *tptr;
                        /* Bump passed bottom pixel because we know it will not */
                        /* be a hole.                                           */
                        tptr+=iw2;
                        mptr+=iw2;
                        bptr+=iw2;
                        /* We bump iy once here and then the FOR bumps it again. */
                        iy++;
                    }
                    else{
                        /* Otherwise, bump to the next pixel below. */
                        tptr+=iw;
                        mptr+=iw;
                        bptr+=iw;
                    }
                }
                /* Bump to start of next column. */
                sptr++;
            }
        }

        void GrayToBin(const int thresh, const int less_pix, const int greater_pix,
                    unsigned char *bdata, const int iw, const int ih)
        {
            for (int i = 0; i < iw * ih; ++i) {
                if (bdata[i] >= thresh) {
                    bdata[i] = static_cast<unsigned char>(greater_pix);
                }
                else {
                    bdata[i] = static_cast<unsigned char>(less_pix);
                }
            }
        }

        void PixelizeMap(std::vector<int>& omap, const int iw, const int ih, int *imap, const int mw, const int mh, const int blocksize)
        {
            std::vector<int> blkoffs;
            int bw, bh;
            BlockOffsets(blkoffs, bw, bh, iw, ih, 0, blocksize);

            if ((bw != mw) || (bh != mh)) {
                throw std::invalid_argument("PixelizeMap: block dimensions do not match");
            }

            omap.resize(static_cast<std::size_t>(iw * ih));
            int* pmap = omap.data();
            for (int bi = 0; bi < mw*mh; ++bi) {
                int* spptr = pmap + blkoffs[static_cast<std::size_t>(bi)];
                for (int y = 0; y < blocksize; ++y) {
                    int* pptr = spptr;
                    for (int x = 0; x < blocksize; ++x) {
                        *pptr++ = imap[static_cast<std::size_t>(bi)];
                    }
                    spptr += iw;
                }
            }
        }

        bool FreePath(const int x1, const int y1, const int x2, const int y2,
                    unsigned char *bdata, const int iw,
                    const LQMParams& lqmParams)
        {
            /* Compute points along line segment between the two points. */
            std::vector<std::pair<int, int>> points;
            LinePoints(points, x1, y1, x2, y2);

            /* Intialize the number of transitions to 0. */
            int trans = 0;
            /* Get the pixel value of first point along line segment. */
            int preval = *(bdata+(y1*iw)+x1);

            /* Foreach remaining point along line segment ... */
            for (int i = 1; i < static_cast<int>(points.size()); ++i) {
                /* Get pixel value of next point along line segment. */
                int nextval = *(bdata+(points[static_cast<std::size_t>(i)].second*iw) + points[static_cast<std::size_t>(i)].first);

                /* If next pixel value different from previous pixel value ... */
                if (nextval != preval) {
                    /* Then we have detected a transition, so bump counter. */
                    ++trans;
                    /* If number of transitions seen > than threshold (ex. 2) ... */
                    if (trans > lqmParams.maxtrans) {
                        /* Return free path to be false. */
                        return false;
                    }
                    /* Otherwise, maximum number of transitions not yet exceeded. */
                    /* Assign the next pixel value to the previous pixel value.   */
                    preval = nextval;
                }
                /* Otherwise, no transition detected this interation. */
            }

            /* If we get here we did not exceed the maximum allowable number of transitions. */
            /* Return free path to be TRUE. */
            return true;
        }

        void LinePoints(std::vector<std::pair<int, int>>& points, const int x1, const int y1, const int x2, const int y2) {
            /* Compute maximum number of points needed to hold line segment. */
            int asize = std::max(std::abs(x2-x1)+2, std::abs(y2-y1)+2);

            /* Allocate point lists to length 'asize'. */
            points.reserve(static_cast<std::size_t>(asize));

            /* Compute delta x and y. */
            int dx = x2 - x1;
            int dy = y2 - y1;

            /* Set x and y increments. */
            int x_incr = (dx >= 0 ? 1 : -1);
            int y_incr = (dy >= 0 ? 1 : -1);

            /* Compute |DX| and |DY|. */
            int adx = std::abs(dx);
            int ady = std::abs(dy);

            /* Set x and y orientations. */
            int inx = adx > ady;
            int iny = ady > adx;

            /*  CASE 1: |DX| > |DY|              */
            /*     Increment in X by +-1         */
            /*               in Y by +-|DY|/|DX| */
            /*        inx   =  1                 */
            /*        iny   =  0                 */
            /*        intx  =  1 (inx)           */
            /*        inty  =  0 (iny)           */
            /*  CASE 2: |DX| < |DY|              */
            /*     Increment in Y by +-1         */
            /*               in X by +-|DX|/|DY| */
            /*        inx   =  0                 */
            /*        iny   =  1                 */
            /*        intx  =  0 (inx)           */
            /*        inty  =  1 (iny)           */
            /*  CASE 3: |DX| == |DY|             */
            /*        inx   =  0                 */
            /*        iny   =  0                 */
            /*        intx  =  1                 */
            /*        inty  =  1                 */
            int intx = 1 - iny;
            int inty = 1 - inx;

            /*                                        DX           */
            /* x_factor = (inx * +-1) +  ( iny * ------------ )    */
            /*                                   max(1, |DY|)      */
            /*                                                     */
            double x_factor = (inx * x_incr) + (iny * (static_cast<double>(dx) / std::max(1, ady)));

            /*                                        DY           */
            /* y_factor = (iny * +-1) +  ( inx * ------------ )    */
            /*                                   max(1, |DX|)      */
            /*                                                     */
            double y_factor = (iny * y_incr) + (inx * (static_cast<double>(dy) / std::max(1, adx)));

            /* Initialize integer coordinates. */
            int ix = x1;
            int iy = y1;
            /* Set floating point coordinates. */
            double rx = static_cast<double>(x1);
            double ry = static_cast<double>(y1);

            /* Assign first point into coordinate list. */
            points.emplace_back(x1, y1);

            while (ix != x2 || iy != y2) {

                if (static_cast<int>(points.size()) >= asize) {
                    throw std::runtime_error("LinePoints generated more points than expected");
                }

                rx += x_factor;
                ry += y_factor;

                /* Need to truncate precision so that answers are consistent */
                /* on different computer architectures when truncating doubles. */
                rx = trunc_dbl_precision(rx);
                ry = trunc_dbl_precision(ry);

                /* Compute new x and y-pixel coords in floating point and  */
                /* then round to the nearest integer.                      */
                ix = (intx * (ix + x_incr)) + (iny * static_cast<int>(rx + 0.5));
                iy = (inty * (iy + y_incr)) + (inx * static_cast<int>(ry + 0.5));

                /* Assign first point into coordinate list. */
                points.emplace_back(ix, iy);
            }
        }

        OpenLQM::Fingerprint CreateNormFingerprint(const OpenLQM::Fingerprint& img) {
            unsigned int ppi = static_cast<unsigned int>(img.resolution);
            unsigned int resFactor = ppi / 500U;
            unsigned int newWidth = WordAlign(static_cast<unsigned int>(std::lrint(static_cast<double>(img.width) / static_cast<double>(resFactor))), false);
            unsigned int newHeight = WordAlign(static_cast<unsigned int>(std::lrint(static_cast<double>(img.height) / static_cast<double>(resFactor))), false);
            unsigned int newSize = newWidth * newHeight;

            OpenLQM::Fingerprint normBB;
            normBB.width = newWidth;
            normBB.height = newHeight;
            normBB.resolution = OpenLQM::PixelDensity::ppi500;
            /* FIXME: https://github.com/usnistgov/openlqm/issues/7 */
            if (img.bitDepth != PixelBitDepth::depth8 || img.bitsPerPixel != 8)
                throw std::runtime_error{"Unsupported bit depth"};
            normBB.bitDepth = img.bitDepth;
            normBB.bitsPerPixel = img.bitsPerPixel;

            normBB.buffer.resize(newSize);
            if (resFactor == 1) {
                memcpy(normBB.buffer.data(), img.buffer.data(), newSize);
            } else {
                for (unsigned int y = 0; y < normBB.height; ++y) {
                    for (unsigned int x = 0; x < normBB.width; ++x) {
                        normBB.buffer[y * newWidth + x] = img.buffer[(y * resFactor * img.width) + (x * resFactor)];
                    }
                }
            }

            return normBB;
        }

        float ClampResolution(float ppi) {
            if (ppi >= 490.0f && ppi <= 510.0f) {
                return 500.0f;
            } else if (ppi >= 980.0f && ppi <= 1020.0f) {
                return 1000.0f;
            } else if (ppi >= 1960.0f && ppi <= 2040.0f) {
                return 2000.0f;
            } else if (ppi >= 3920.0f && ppi <= 4080.0f) {
                return 4000.0f;
            }
            return ppi;
        }
    }
}
