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

#include <openlqm/openlqm_minutia.hpp>
#include <openlqm/openlqm_img_util.hpp>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <sstream>
#include <array>
#include <cstring>

namespace OpenLQM {
    namespace Core {
        Shape::Shape(const int *contour_x, const int *contour_y, const int ncontour)
        {
            /* Find xmin, ymin, xmax, ymax on contour. */

            auto xMinMax = std::minmax_element(contour_x, contour_x + ncontour);
            auto yMinMax = std::minmax_element(contour_y, contour_y + ncontour);

            int xmin = *xMinMax.first;
            int ymin_ = *yMinMax.first;
            int xmax = *xMinMax.second;
            int ymax_ = *yMinMax.second;

            /* Allocate and initialize a shape structure. */

            /* Compute allocation parameters. */
            /* First, compute the number of scanlines spanned by the shape. */
            int alloc_rows = ymax_ - ymin_ + 1;
            /* Second, compute the "maximum" number of contour points possible    */
            /* on a row.  Here we are allocating the maximum number of contiguous */
            /* pixels on each row which will be sufficiently larger than the      */
            /* number of actual contour points.                                   */
            int alloc_pts = xmax - xmin + 1;

            /* Allocate the list of rows.  We now this number will fit */
            /* the shape exactly.                                      */

            this->ymin = ymin_;
            this->ymax = ymax_;

            this->rows.reserve(static_cast<std::size_t>(alloc_rows));
            for (int i = 0, y = ymin_; i < alloc_rows; ++i, ++y) {
                this->rows.emplace_back(y, alloc_pts);

            }

            /* Foreach point on contour ... */
            for (int i = 0; i < ncontour; ++i) {
                /* Add point to corresponding row. */
                /* First set a pointer to the current row.  We need to subtract */
                /* ymin_ because the rows are indexed relative to the top-most   */
                /* scanline in the shape.                                       */
                Row& row = this->rows[static_cast<std::size_t>(contour_y[static_cast<std::size_t>(i)] - ymin_)];

                /* It is possible with complex shapes to reencounter points        */
                /* already visited on a contour, especially at "pinching" points   */
                /* along the contour.  So we need to test to see if a point has    */
                /* already been stored in the row.  If not in row list already ... */
                if (std::find(row.xs.begin(), row.xs.end(), contour_x[static_cast<std::size_t>(i)]) == row.xs.end()) {
                    /* Assign the x-coord of the current contour point to the row */
                    /* and bump the row's point counter.  All the contour points  */
                    /* on the same row share the same y-coord.                    */
                    row.xs.push_back(contour_x[static_cast<std::size_t>(i)]);
                }
                /* Otherwise, point is already stored in row, so ignore. */
            }

            /* Foreach row in the shape. */
            for (Row& row : this->rows) {
                /* Sort row points increasing on their x-coord. */
                std::sort(row.xs.begin(), row.xs.end());
            }
        }

        Dir2Rad::Dir2Rad(int ndirs_) : ndirs(ndirs_), cosLut(static_cast<std::size_t>(ndirs_)), sinLut(static_cast<std::size_t>(ndirs_)) {
            /* Pi_factor sets the period of the trig functions to NDIRS units in x. */
            /* For example, if NDIRS==16, then pi_factor = 2(PI/16) = .3926...      */
            double pi_factor = 2.0*LQM_PI/static_cast<double>(ndirs);

            /* Now compute cos and sin values for each direction.    */
            for (int i = 0; i < ndirs; ++i) {
                double theta = static_cast<double>(i * pi_factor);
                double cs = std::cos(theta);
                double sn = std::sin(theta);
                /* Need to truncate precision so that answers are consistent */
                /* on different computer architectures. */
                cs = trunc_dbl_precision(cs);
                sn = trunc_dbl_precision(sn);
                cosLut[static_cast<std::size_t>(i)] = cs;
                sinLut[static_cast<std::size_t>(i)] = sn;
            }
        }

        void DFTWave::Init(int blocksize, double freq) {
            cosLut.resize(static_cast<std::size_t>(blocksize));
            sinLut.resize(static_cast<std::size_t>(blocksize));

            for (int j = 0; j < blocksize; j++) {
                /* Compute sample points from frequency */
                double x = freq * static_cast<double>(j);
                /* Store cos and sin components of sample point */
                cosLut[static_cast<std::size_t>(j)] = std::cos(x);
                sinLut[static_cast<std::size_t>(j)] = std::sin(x);
            }
        }

        DFTWaves::DFTWaves(int nwaves_, int blocksize) : nwaves(nwaves_), wavelen(blocksize), waves(static_cast<std::size_t>(nwaves_)) {
            /* Pi_factor sets the period of the trig functions to BLOCKSIZE units */
            /* in x.  For example, if BLOCKSIZE==24, then                         */
            /*                         pi_factor = 2(PI/24) = .26179...           */
            double pi_factor = 2.0*LQM_PI/static_cast<double>(blocksize);

            /* Foreach of 4 DFT frequency coef ... */
            for (int i = 0; i < nwaves; ++i) {
                /* Compute actual frequency */
                double freq = pi_factor * OpenLQM::DFT_COEFFICIENTS[i];
                waves[static_cast<std::size_t>(i)].Init(blocksize, freq);
            }
        }

        RotGrids::RotGrids(
            const int iw, const int ipad,
            const double start_dir_angle, const int ndirs,
            const int grid_w_, const int grid_h_, const GridOffsetMode relative2_
        ) : relative2(relative2_), start_angle(start_dir_angle), ngrids(ndirs), grid_w(grid_w_), grid_h(grid_h_)
        {
            double diag = sqrt(static_cast<double>((grid_w*grid_w)+(grid_h*grid_h)));
            int grid_pad, min_dim;

            switch (relative2) {
                case RELATIVE2CENTER:
                    /* Assumption: all grid centers reside in valid/allocated memory. */
                    this->pad = static_cast<int>((diag-1)/static_cast<double>(2.0));
                    /* Need to truncate precision so that answers are consistent */
                    /* on different computer architectures when rounding doubles. */
                    this->pad = static_cast<int>(trunc_dbl_precision(this->pad));
                    grid_pad = sround(this->pad);
                    break;
                case RELATIVE2ORIGIN:
                    /* Assumption: all grid origins reside in valid/allocated memory. */
                    min_dim = std::min(grid_w, grid_h);
                    /* Compute pad as difference between the smallest grid dimension */
                    /* and the diagonal distance of the grid. */
                    this->pad = static_cast<int>((diag - min_dim) / static_cast<double>(2.0));
                    /* Need to truncate precision so that answers are consistent */
                    /* on different computer architectures when rounding doubles. */
                    this->pad = static_cast<int>(trunc_dbl_precision(this->pad));
                    grid_pad = sround(this->pad);
                    break;
                default:
                    throw std::invalid_argument(std::string("RotGrids: Illegal relative flag : ") + std::to_string(static_cast<int>(relative2)));
            }

            if (ipad == LQM_UNDEFINED) {
                /* Use the padding specifically required by the rotated grids herein. */
                this->pad = grid_pad;
            } else {
                /* Otherwise, input pad was specified, so check to make sure it is */
                /* sufficiently large to handle the rotated grids herein.          */
                if (ipad < grid_pad) {
                    /* ERROR: input pad is NOT large enough */
                    throw std::invalid_argument("RotGrids: Pad passed is too small");
                }
                /* Otherwise, use the specified input pad in computing grid offsets. */
                this->pad = ipad;
            }

            /* Total number of points in grid */
            int grid_size = grid_w * grid_h;

            /* Compute width of "padded" image */
            int pw = iw + (this->pad<<1);

            /* Center coord of grid (0-oriented). */
            double cx = (grid_w-1)/static_cast<double>(2.0);
            double cy = (grid_h-1)/static_cast<double>(2.0);

            /* Allocate list of rotgrid pointers */
            this->grids.resize(static_cast<std::size_t>(ndirs));

            /* Pi_offset is the offset in radians from which angles are to begin. */
            double pi_offset = start_dir_angle;
            double pi_incr = LQM_PI/static_cast<double>(ndirs);  /* if ndirs == 16, incr = 11.25 degrees */

            double theta = pi_offset;
            for (int dir = 0; dir < ndirs; dir++, theta += pi_incr) {
                /* Allocate a rotgrid */
                this->grids[static_cast<std::size_t>(dir)].resize(static_cast<std::size_t>(grid_size));
                int* pGrid = this->grids[static_cast<std::size_t>(dir)].data();

                /* Compute cos and sin of current angle */
                double cs = std::cos(theta);
                double sn = std::sin(theta);

                /* This next section of nested FOR loops precomputes a         */
                /* rotated grid.  The rotation is set up to rotate a GRID_W X  */
                /* GRID_H grid on its center point at C=(Cx,Cy). The current   */
                /* pixel being rotated is P=(Ix,Iy).  Therefore, we have a     */
                /* rotation transformation of point P about pivot point C.     */
                /* The rotation transformation about a pivot point in matrix   */
                /* form is:                                                    */
                /*
                        +-                                                       -+
                        |             cos(T)                   sin(T)           0 |
            [Ix Iy 1] |            -sin(T)                   cos(T)           0 |
                        | (1-cos(T))*Cx + Cy*sin(T)  (1-cos(T))*Cy - Cx*sin(T)  1 |
                        +-                                                       -+
                */
                /* Multiplying the 2 matrices and combining terms yeilds the */
                /* equations for rotated coordinates (Rx, Ry):               */
                /*        Rx = Cx + (Ix - Cx)*cos(T) - (Iy - Cy)*sin(T)      */
                /*        Ry = Cy + (Ix - Cx)*sin(T) + (Iy - Cy)*cos(T)      */
                /*                                                           */
                /* Care has been taken to ensure that (for example) when     */
                /* BLOCKSIZE==24 the rotated indices stay within a centered  */
                /* 34X34 area.                                               */
                /* This is important for computing an accurate padding of    */
                /* the input image.  The rotation occurs "in-place" so that  */
                /* outer pixels in the grid are mapped at times from         */
                /* adjoining blocks.  As a result, to keep from accessing    */
                /* "unknown" memory or pixels wrapped from the other side of */
                /* the image, the input image should first be padded by      */
                /* PAD=round((DIAG - BLOCKSIZE)/2.0) where DIAG is the       */
                /* diagonal distance of the grid.                            */
                /* For example, when BLOCKSIZE==24, Dx=34, so PAD=5.         */

                /* Foreach each y coord in block ... */
                for (int iy = 0; iy < grid_h; ++iy) {
                    /* Compute rotation factors dependent on Iy (include constant) */
                    double fxm = -1.0 * ((iy - cy) * sn);
                    double fym = ((iy - cy) * cs);

                    /* If offsets are to be relative to the grids origin, then */
                    /* we need to subtract CX and CY.                          */
                    if (relative2 == RELATIVE2ORIGIN) {
                        fxm += cx;
                        fym += cy;
                    }

                    /* foreach each x coord in block ... */
                    for (int ix = 0; ix < grid_w; ++ix) {

                        /* Now combine factors dependent on Iy with those of Ix */
                        double fx = fxm + ((ix - cx) * cs);
                        double fy = fym + ((ix - cx) * sn);
                        /* Need to truncate precision so that answers are consistent */
                        /* on different computer architectures when rounding doubles. */
                        fx = trunc_dbl_precision(fx);
                        fy = trunc_dbl_precision(fy);
                        int ixt = sround(fx);
                        int iyt = sround(fy);

                        /* Store the current pixels relative   */
                        /* rotated offset.  Make sure to       */
                        /* multiply the y-component of the     */
                        /* offset by the "padded" image width! */
                        *pGrid++ = ixt + (iyt * pw);
                    } /* ix */
                } /* iy */
            } /* dir */
        }

        void DetectMinutiae(
            Minutiae& ominutiae,
            std::vector<int>& direction_map, std::vector<int>& low_contrast_map, std::vector<int>& low_flow_map, std::vector<int>& high_curve_map,
            int *omw, int *omh,
            std::vector<unsigned char>& binarizedData, int *obw, int *obh, // median grayscale
            unsigned char *idata, const int iw, const int ih,
            const LQMParams& lqmParams,
            unsigned char *gray10pctMap,  // grayscale of 10th %ile (dark)
            unsigned char *gray90pctMap,  // grayscale of 90th %ile (light)
            unsigned char *grayMedianMap,  // median grayscale
            unsigned char *grayRangeMap,   // 90%-10% gray value
            unsigned char *grayCountMap,   // number of gray values in use
            double *maxMagnitude,
            double *normMagnitude,
            double *lowFreqMagnitude,
            unsigned char *directionChangeMap,// changes in direction map from initial to final versions
            unsigned char *validNeighbors,
            unsigned char *curvatureMap
        ) {
            int maxpad = OpenLQM::Core::GetMaxPadding(lqmParams.windowsize, lqmParams.windowoffset, lqmParams.dirbin_grid_w, lqmParams.dirbin_grid_h);

            /* Initialize lookup table for converting integer directions */
            Dir2Rad dir2rad(lqmParams.num_directions);

            /* Initialize wave 4form lookup tables for DFT analyses. */
            /* used for direction binarization.                             */
            DFTWaves dftwaves(lqmParams.num_dft_waves, lqmParams.windowsize);

            /* Initialize lookup table for pixel offsets to rotated grids */
            /* used for DFT analyses.                                     */
            RotGrids dftgrids(iw, maxpad,
                lqmParams.start_dir_angle, lqmParams.num_directions,
                lqmParams.windowsize, lqmParams.windowsize,
                GridOffsetMode::RELATIVE2ORIGIN
            );

            std::vector<unsigned char> paddedData;

            /* Pad input image based on max padding. */
            int pw, ph;
            if (maxpad > 0) {   /* May not need to pad at all */
                OpenLQM::Core::PadUCharImage(paddedData, pw, ph, idata, iw, ih,
                                    maxpad, lqmParams.pad_value);
            }
            else {
                /* If padding is unnecessary, then copy the input image. */
                paddedData.resize(static_cast<std::size_t>(iw * ih));
                memcpy(paddedData.data(), idata, static_cast<std::size_t>(iw * ih));
                pw = iw;
                ph = ih;
            }

            /* Scale input image to 6 bits [0..63] */
            /* !!! Would like to remove this dependency eventualy !!!     */
            /* But, the DFT computations will need to be changed, and     */
            /* could not get this work upon first attempt. Also, if not   */
            /* careful, I think accumulated power magnitudes may overflow */
            /* doubles.                                                   */

            OpenLQM::Core::Bits8to6(paddedData.data(), pw, ph);

            /******************/
            /*      MAPS      */
            /******************/

            /* Generate block maps from the input image. */
            int mw, mh;
            OpenLQM::Core::GenImageMaps(direction_map, low_contrast_map,
                                low_flow_map, high_curve_map, &mw, &mh,
                                paddedData.data(), pw, ph, dir2rad, dftwaves, dftgrids, lqmParams,
                                gray10pctMap, gray90pctMap, grayMedianMap, grayRangeMap, grayCountMap,
                                maxMagnitude, normMagnitude, lowFreqMagnitude,
                                directionChangeMap, validNeighbors, curvatureMap);

            /******************/
            /* BINARIZATION   */
            /******************/

            /* Initialize lookup table for pixel offsets to rotated grids */
            /* used for directional binarization.                         */
            RotGrids dirbingrids(
                iw, maxpad,
                lqmParams.start_dir_angle, lqmParams.num_directions,
                lqmParams.dirbin_grid_w, lqmParams.dirbin_grid_h,
                GridOffsetMode::RELATIVE2CENTER
            );

            /* Binarize input image based on NMAP information. */
            int bw, bh;
            OpenLQM::Core::Binarize(
                binarizedData, bw, bh,
                paddedData.data(), pw, ph, direction_map.data(), mw,
                dirbingrids, lqmParams
            );

            /* Check dimension of binary image.  If they are different from */
            /* the input image, then raise exception                        */
            if (iw != bw || ih != bh) {
                std::stringstream ss;
                ss << "DetectMinutiae() Invalid argument: binary image has bad dimensions: " << bw << ", " << bh << std::endl;
                throw std::invalid_argument(ss.str());
            }

            /******************/
            /*   DETECTION    */
            /******************/

            /* Convert 8-bit grayscale binary image [0,255] to */
            /* 8-bit binary image [0,1].                       */
            OpenLQM::Core::GrayToBin(1, 1, 0, binarizedData.data(), iw, ih);

            /* Reserve capacity for detected minutiae. */
            ominutiae.reserve(LQM_INITIAL_MINUTIA_CAPACITY);

            /* Detect the minutiae in the binarized image. */
            FindMinutiae(ominutiae, binarizedData.data(), iw, ih,
                                        direction_map.data(), low_flow_map.data(), high_curve_map.data(),
                                        mw, mh, lqmParams);

            RemoveFalseMinutia(ominutiae, binarizedData.data(), iw, ih,
                            direction_map.data(), low_flow_map.data(), high_curve_map.data(), mw, mh,
                            lqmParams);

            /******************/
            /*  RIDGE COUNTS  */
            /******************/

            CountMinutiaeRidges(ominutiae, binarizedData.data(), iw, ih, lqmParams);

            /******************/
            /*    WRAP-UP     */
            /******************/

            /* Convert 8-bit binary image [0,1] to 8-bit */
            /* grayscale binary image [0,255].           */
            OpenLQM::Core::GrayToBin(1, 255, 0, binarizedData.data(), iw, ih);

            *omw = mw;
            *omh = mh;
            *obw = bw;
            *obh = bh;
        }

        void FindMinutiae(Minutiae& minutiae,
                    unsigned char *bdata, const int iw, const int ih,
                    int *direction_map, int *low_flow_map, int *high_curve_map,
                    const int mw, const int mh,
                    const LQMParams& lqmParams)
        {
            /* Pixelize the maps by assigning block values to individual pixels. */

            std::vector<int> pixDirectionMap;
            OpenLQM::Core::PixelizeMap(pixDirectionMap, iw, ih, direction_map, mw, mh, lqmParams.blocksize);
            std::vector<int> pixLowFlowMap;
            OpenLQM::Core::PixelizeMap(pixLowFlowMap, iw, ih, low_flow_map, mw, mh, lqmParams.blocksize);
            std::vector<int> pixHighCurveMap;
            OpenLQM::Core::PixelizeMap(pixHighCurveMap, iw, ih, high_curve_map, mw, mh, lqmParams.blocksize);

            ScanForMinutiaeHorizontally(minutiae, bdata, iw, ih, pixDirectionMap.data(), pixLowFlowMap.data(), pixHighCurveMap.data(), lqmParams);
            ScanForMinutiaeVertically(minutiae, bdata, iw, ih, pixDirectionMap.data(), pixLowFlowMap.data(), pixHighCurveMap.data(), lqmParams);
        }

        void ScanForMinutiaeHorizontally(Minutiae& minutiae,
                        unsigned char *bdata, const int iw, const int ih,
                        int *pdirection_map, int *plow_flow_map, int *phigh_curve_map,
                        const LQMParams& lqmParams)
        {
            int possible[LQM_NFEATURES], nposs;

            /* Set scan region to entire image. */
            int sx = 0;
            int ex = iw;
            int sy = 0;
            int ey = ih;

            /* Start at first row in region. */
            int cy = sy;
            /* While second scan row not outside the bottom of the scan region... */
            while (cy+1 < ey) {
                /* Start at beginning of new scan row in region. */
                int cx = sx;
                /* While not at end of region's current scan row. */
                while (cx < ex) {
                    /* Get pixel pair from current x position in current and next */
                    /* scan rows. */
                    unsigned char* p1ptr = bdata+(cy*iw)+cx;
                    unsigned char* p2ptr = bdata+((cy+1)*iw)+cx;
                    /* If scan pixel pair matches first pixel pair of */
                    /* 1 or more features... */
                    if (MatchFirstPair(*p1ptr, *p2ptr, possible, &nposs)) {
                        /* Bump forward to next scan pixel pair. */
                        cx++;
                        p1ptr++;
                        p2ptr++;
                        /* If not at end of region's current scan row... */
                        if(cx < ex){
                        /* If scan pixel pair matches second pixel pair of */
                        /* 1 or more features... */
                        if (MatchSecondPair(*p1ptr, *p2ptr, possible, &nposs)) {
                            /* Store current x location. */
                            int x2 = cx;
                            /* Skip repeated pixel pairs. */
                            SkipRepeatedHorizontalPair(&cx, ex, &p1ptr, &p2ptr);
                            /* If not at end of region's current scan row... */
                            if (cx < ex) {
                                /* If scan pixel pair matches third pixel pair of */
                                /* a single feature... */
                                if (MatchThirdPair(*p1ptr, *p2ptr, possible, &nposs)) {
                                    /* Process detected minutia point. */
                                    ProcessHorizontalScanMinutia(minutiae,
                                        cx, cy, x2, possible[0],
                                        bdata, iw, ih, pdirection_map,
                                        plow_flow_map, phigh_curve_map,
                                        lqmParams);
                                    /* Return code may be:                         */
                                    /* 1.  ret==0 (normal)                         */
                                    /* 2. ret==LQM_IGNORE (ignore current feature) */
                                }

                                /* Set up to resume scan. */
                                /* Test to see if 3rd pair can slide into 2nd pair. */
                                /* The values of the 2nd pair MUST be different.    */
                                /* If 3rd pair values are different ... */
                                if (*p1ptr != *p2ptr) {
                                    /* Set next first pair to last of repeated */
                                    /* 2nd pairs, ie. back up one pair.        */
                                    cx--;
                                }

                                /* Otherwise, 3rd pair can't be a 2nd pair, so  */
                                /* keep pointing to 3rd pair so that it is used */
                                /* in the next first pair test.                 */

                            } /* Else, at end of current scan row. */
                        }

                        /* Otherwise, 2nd pair failed, so keep pointing to it */
                        /* so that it is used in the next first pair test.    */

                        } /* Else, at end of current scan row. */
                    }
                    /* Otherwise, 1st pair failed... */
                    else{
                        /* Bump forward to next pixel pair. */
                        cx++;
                    }
                } /* While not at end of current scan row. */
                /* Bump forward to next scan row. */
                cy++;
            } /* While not out of scan rows. */
        }

        int MatchFirstPair(unsigned char p1, unsigned char p2, int *possible, int *nposs)
        {
            /* Set possibilities to 0 */
            *nposs = 0;

            /* Foreach set of feature pairs ... */
            for (int i = 0; i < LQM_NFEATURES; ++i) {
                /* If current scan pair matches first pair for feature ... */
                if (p1 == OpenLQM::FEATURE_PATTERNS[i].first[0] && p2 == OpenLQM::FEATURE_PATTERNS[i].first[1]) {
                    /* Store feature as a possible match. */
                    possible[*nposs] = i;
                    /* Bump number of stored possibilities. */
                    (*nposs)++;
                }
            }

            /* Return number of stored possibilities. */
            return *nposs;
        }

        int MatchSecondPair(unsigned char p1, unsigned char p2, int *possible, int *nposs)
        {
            /* Store input possibilities. */
            int tnposs = *nposs;
            /* Reset output possibilities to 0. */
            *nposs = 0;

            /* If current scan pair values are the same ... */
            if (p1 == p2) {
                /* Simply return because pair can't be a second feature pair. */
                return *nposs;
            }

            /* Foreach possible match based on first pair ... */
            for (int i = 0; i < tnposs; i++) {
                /* If current scan pair matches second pair for feature ... */
                if (p1 == OpenLQM::FEATURE_PATTERNS[possible[i]].second[0] && p2 == OpenLQM::FEATURE_PATTERNS[possible[i]].second[1]) {
                    /* Store feature as a possible match. */
                    possible[*nposs] = possible[i];
                    /* Bump number of stored possibilities. */
                    (*nposs)++;
                }
            }

            /* Return number of stored possibilities. */
            return *nposs;
        }

        void SkipRepeatedHorizontalPair(int *cx, const int ex,
                        unsigned char **p1ptr, unsigned char **p2ptr)
        {
            /* Store starting pixel pair. */
            int old1 = **p1ptr;
            int old2 = **p2ptr;

            /* Bump horizontally to next pixel pair. */
            (*cx)++;
            (*p1ptr)++;
            (*p2ptr)++;

            /* While not at right of scan region... */
            while (*cx < ex) {
                /* If one or the other pixels in the new pair are different */
                /* from the starting pixel pair...                          */
                if ((**p1ptr != old1) || (**p2ptr != old2)) {
                    /* Done skipping repreated pixel pairs. */
                    return;
                }
                /* Otherwise, bump horizontally to next pixel pair. */
                (*cx)++;
                (*p1ptr)++;
                (*p2ptr)++;
            }
        }

        void SkipRepeatedVerticalPair(int *cy, const int ey,
                        unsigned char **p1ptr, unsigned char **p2ptr,
                        const int iw)
        {
            /* Store starting pixel pair. */
            int old1 = **p1ptr;
            int old2 = **p2ptr;

            /* Bump vertically to next pixel pair. */
            (*cy)++;
            (*p1ptr)+=iw;
            (*p2ptr)+=iw;

            /* While not at bottom of scan region... */
            while(*cy < ey) {
                /* If one or the other pixels in the new pair are different */
                /* from the starting pixel pair...                          */
                if ((**p1ptr != old1) || (**p2ptr != old2)) {
                    /* Done skipping repreated pixel pairs. */
                    return;
                }
                /* Otherwise, bump vertically to next pixel pair. */
                (*cy)++;
                (*p1ptr)+=iw;
                (*p2ptr)+=iw;
            }
        }

        int MatchThirdPair(unsigned char p1, unsigned char p2, int *possible, int *nposs)
        {
            /* Store input possibilities. */
            int tnposs = *nposs;
            /* Reset output possibilities to 0. */
            *nposs = 0;

            /* Foreach possible match based on first and second pairs ... */
            for (int i = 0; i < tnposs; i++) {
                /* If current scan pair matches third pair for feature ... */
                if (p1 == OpenLQM::FEATURE_PATTERNS[possible[i]].third[0] && p2 == OpenLQM::FEATURE_PATTERNS[possible[i]].third[1]) {
                    /* Store feature as a possible match. */
                    possible[*nposs] = possible[i];
                    /* Bump number of stored possibilities. */
                    (*nposs)++;
                }
            }

            /* Return number of stored possibilities. */
            return *nposs;
        }

        int ProcessHorizontalScanMinutia(Minutiae& minutiae,
                        const int cx, const int cy,
                        const int x2, const int feature_id,
                        unsigned char *bdata, const int iw, const int ih,
                        int *pdirection_map, int *plow_flow_map, int *phigh_curve_map,
                        const LQMParams& lqmParams)
        {
            int idir;
            int dmapval, fmapval, cmapval;
            double reliability;

            /* Set x location of minutia point to be half way between */
            /* first position of second feature pair and position of  */
            /* third feature pair.                                    */
            int x_loc = (cx + x2)>>1;
            int y_loc;

            /* Set same x location to neighboring edge pixel. */
            int x_edge = x_loc;
            int y_edge;

            /* Feature location should always point to either ending  */
            /* of ridge or (for bifurcations) ending of valley.       */
            /* So, if detected feature is APPEARING...                */
            if (OpenLQM::FEATURE_PATTERNS[feature_id].appearing) {
                /* Set y location to second scan row. */
                y_loc = cy+1;
                /* Set y location of neighboring edge pixel to the first scan row. */
                y_edge = cy;
            }
            /* Otherwise, feature is DISAPPEARING... */
            else {
                /* Set y location to first scan row. */
                y_loc = cy;
                /* Set y location of neighboring edge pixel to the second scan row. */
                y_edge = cy+1;
            }

            dmapval = *(pdirection_map+(y_loc*iw)+x_loc);
            fmapval = *(plow_flow_map+(y_loc*iw)+x_loc);
            cmapval = *(phigh_curve_map+(y_loc*iw)+x_loc);

            /* If the minutia point is in a block with INVALID direction ... */
            if (dmapval == LQM_INVALID_DIR) {
                /* Then, IGNORE the point. */
                return LQM_IGNORE;
            }

            /* If current minutia is in a HIGH CURVATURE block ... */
            if (cmapval) {
                /* Adjust location and direction locally. */
                int adjustRet = AdjustHighCurvatureMinutia(&idir, &x_loc, &y_loc,
                                    &x_edge, &y_edge, x_loc, y_loc, x_edge, y_edge,
                                    bdata, iw, ih, plow_flow_map, minutiae, lqmParams);
                if (adjustRet) {
                    /* Returned LQM_IGNORE, so forward */
                    return adjustRet;
                }
                /* Otherwise, we have our high-curvature minutia attributes. */
            }
            /* Otherwise, minutia is in fairly low-curvature block... */
            else {
                /* Get minutia direction based on current block's direction. */
                idir = GetLowCurvatureDirection(LQM_SCAN_HORIZONTAL,
                                OpenLQM::FEATURE_PATTERNS[feature_id].appearing, dmapval,
                                lqmParams.num_directions);
            }

            /* If current minutia is in a LOW RIDGE FLOW block ... */
            if (fmapval) {
                reliability = OpenLQM::Presets::MinutiaeReliability::MEDIUM;
            }
            else {
                /* Otherwise, minutia is in a block with reliable direction and */
                /* binarization.                                                */
                reliability = OpenLQM::Presets::MinutiaeReliability::HIGH;
            }

            {
                /* Create a minutia object based on derived attributes. */
                Minutia minutia(
                    x_loc, y_loc, x_edge, y_edge, idir,
                    reliability,
                    OpenLQM::FEATURE_PATTERNS[feature_id].type,
                    OpenLQM::FEATURE_PATTERNS[feature_id].appearing, feature_id
                );

                /* Update the minutiae list with potential new minutia. */
                ScanAndUpdateMinutiae(minutiae, minutia, LQM_SCAN_HORIZONTAL,
                                            dmapval, bdata, iw, ih, lqmParams);
            }

            /* Return normally */
            return 0;
        }

        int ProcessVerticalScanMinutia(Minutiae& minutiae,
                        const int cx, const int cy,
                        const int y2, const int feature_id,
                        unsigned char *bdata, const int iw, const int ih,
                        int *pdirection_map, int *plow_flow_map, int *phigh_curve_map,
                        const LQMParams& lqmParams)
        {
            int x_loc, x_edge;
            int idir, ret;

            /* Feature location should always point to either ending  */
            /* of ridge or (for bifurcations) ending of valley.       */
            /* So, if detected feature is APPEARING...                */
            if (OpenLQM::FEATURE_PATTERNS[feature_id].appearing) {
                /* Set x location to second scan column. */
                x_loc = cx+1;
                /* Set x location of neighboring edge pixel to the first scan column. */
                x_edge = cx;
            }
            /* Otherwise, feature is DISAPPEARING... */
            else {
                /* Set x location to first scan column. */
                x_loc = cx;
                /* Set x location of neighboring edge pixel to the second scan column. */
                x_edge = cx+1;
            }

            /* Set y location of minutia point to be half way between */
            /* first position of second feature pair and position of  */
            /* third feature pair.                                    */
            int y_loc = (cy + y2)>>1;
            /* Set same y location to neighboring edge pixel. */
            int y_edge = y_loc;

            int dmapval = *(pdirection_map+(y_loc*iw)+x_loc);
            int fmapval = *(plow_flow_map+(y_loc*iw)+x_loc);
            int cmapval = *(phigh_curve_map+(y_loc*iw)+x_loc);

            /* If the minutia point is in a block with INVALID direction ... */
            if(dmapval == LQM_INVALID_DIR) {
                /* Then, IGNORE the point. */
                return LQM_IGNORE;
            }

            /* If current minutia is in a HIGH CURVATURE block... */
            if(cmapval) {
                /* Adjust location and direction locally. */
                ret = AdjustHighCurvatureMinutia(&idir, &x_loc, &y_loc,
                    &x_edge, &y_edge, x_loc, y_loc, x_edge, y_edge,
                    bdata, iw, ih, plow_flow_map, minutiae, lqmParams);
                if (ret) {
                    /* IGNORE minutia. */
                    return ret;
                }
                /* Otherwise, we have our high-curvature minutia attributes. */
            }
            /* Otherwise, minutia is in fairly low-curvature block... */
            else {
                /* Get minutia direction based on current block's direction. */
                idir = GetLowCurvatureDirection(LQM_SCAN_VERTICAL,
                                OpenLQM::FEATURE_PATTERNS[feature_id].appearing, dmapval,
                                lqmParams.num_directions);
            }

            double reliability;
            /* If current minutia is in a LOW RIDGE FLOW block ... */
            if(fmapval) {
                reliability = OpenLQM::Presets::MinutiaeReliability::MEDIUM;
            }
            else {
                /* Otherwise, minutia is in a block with reliable direction and */
                /* binarization.                                                */
                reliability = OpenLQM::Presets::MinutiaeReliability::HIGH;
            }

            /* Create a minutia object based on derived attributes. */
            {
                Minutia minutia(
                    x_loc, y_loc, x_edge, y_edge, idir,
                    reliability,
                    OpenLQM::FEATURE_PATTERNS[feature_id].type,
                    OpenLQM::FEATURE_PATTERNS[feature_id].appearing, feature_id
                );

                /* Update the minutiae list with potential new minutia. */
                ScanAndUpdateMinutiae(minutiae, minutia, LQM_SCAN_VERTICAL, dmapval, bdata, iw, ih, lqmParams);
            }

            /* Return normally. */
            return 0;
        }

        int AdjustHighCurvatureMinutia(int *oidir, int *ox_loc, int *oy_loc,
                    int *ox_edge, int *oy_edge,
                    const int x_loc, const int y_loc,
                    const int x_edge, const int y_edge,
                    unsigned char *bdata, const int iw, const int ih,
                    int *plow_flow_map, Minutiae& minutiae, const LQMParams& lqmParams)
        {
            /* Set variable from parameter structure. */
            int half_contour = lqmParams.high_curve_half_contour;

            /* Set edge length for computing contour's angle of curvature */
            /* to one quarter of desired pixel length of entire contour.  */
            /* Ex. If half_contour==14, then contour length==29=(2X14)+1  */
            /* and angle_edge==7=(14/2).                                  */
            int angle_edge = half_contour>>1;

            /* Get the pixel value of current feature. */
            int feature_pix = *(bdata + (y_loc * iw) + x_loc);

            /* Extract feature's contour. */
            Contour contour;
            int ret = GetHighCurvatureContour(contour, half_contour,
                                    x_loc, y_loc, x_edge, y_edge, bdata, iw, ih);
            if (ret) {
                /* Returns with:                                                    */
                /*    1. Successful or empty contour == 0                           */
                /*       If contour is empty, then contour lists are not allocated. */
                /*    2. Contour forms loop == LQM_LOOP_FOUND                       */

                /* If the contour forms a loop... */
                if (ret == LQM_LOOP_FOUND) {
                    /* If the order of the contour is clockwise, then the loops's     */
                    /* contour pixels are outside the corresponding edge pixels.  We  */
                    /* definitely do NOT want to fill based on the feature pixel in   */
                    /* this case, because it is OUTSIDE the loop.  For now we will    */
                    /* ignore the loop and the minutia that triggered its tracing.    */
                    /* It is likely that other minutia on the loop will be            */
                    /* detected that create a contour on the "inside" of the loop.    */
                    /* There is another issue here that could be addressed ...        */
                    /* It seems that many/multiple minutia are often detected within  */
                    /* the same loop, which currently requires retracing the loop,    */
                    /* locating minutia on opposite ends of the major axis of the     */
                    /* loop, and then determining that the minutia have already been  */
                    /* entered into the minutiae list upon processing the very first   */
                    /* minutia detected in the loop.  There is a lot of redundant     */
                    /* work being done here!                                          */
                    /* Is_loop_clockwise takes a default value to be returned if the  */
                    /* routine is unable to determine the direction of the contour.   */
                    /* In this case, we want to IGNORE the loop if we can't tell its  */
                    /* direction so that we do not inappropriately fill the loop, so  */
                    /* we are passing the default value TRUE.                         */
                    ret = IsLoopClockwise(contour.xV.data(), contour.yV.data(), contour.ncontour, LQM_TRUE);
                    if (ret) {
                        /* Loop is clockwise, so return LQM_IGNORE. */
                        return LQM_IGNORE;
                    }

                    /* Otherwise, process the clockwise-ordered contour of the loop */
                    /* as it may contain minutia.  If no minutia found, then it is  */
                    /* filled in.                                                   */
                    ProcessMinutiaeLoop(minutiae, contour, bdata, iw, ih, plow_flow_map, lqmParams);
                    /* Either a minutia pair was extracted or the loop was     */
                    /* filled.  Either way we want to ignore the minutia that   */
                    /* started the whole loop processing in the beginning.      */
                    return LQM_IGNORE;
                }

                /* Otherwise not a loop, so get_high_curvature_contour incurred */
                /* a system error.  Return the error code.                      */
                return ret;
            }

            /* If contour is empty ... then contour lists were not allocated, so */
            /* simply return LQM_IGNORE. The contour comes back empty when there */
            /* were not a sufficient number of points found on the contour.      */
            if (contour.ncontour == 0) {
                return LQM_IGNORE;
            }

            /* Otherwise, there are contour points to process. */

            /* Given the contour, determine the point of highest curvature */
            /* (ie. forming the minimum angle between contour walls).      */
            double min_theta;
            int min_i;
            ret = MinContourTheta(&min_i, &min_theta, angle_edge,
                contour.xV.data(), contour.yV.data(), contour.ncontour);
            if (ret) {
                /* Returned LQM_IGNORE. Forward */
                return ret;
            }

            /* If the minimum theta found along the contour is too large... */
            if (min_theta >= lqmParams.max_high_curve_theta) {
                /* Reject the high-curvature minutia, and return LQM_IGNORE. */
                return LQM_IGNORE;
            }

            /* Test to see if interior of curvature is OK.  Compute midpoint   */
            /* between left and right points symmetrically distant (angle_edge */
            /* pixels) from the contour's point of minimum theta.              */
            int mid_x = (contour.xV[static_cast<std::size_t>(min_i-angle_edge)] + contour.xV[static_cast<std::size_t>(min_i+angle_edge)])>>1;
            int mid_y = (contour.yV[static_cast<std::size_t>(min_i-angle_edge)] + contour.yV[static_cast<std::size_t>(min_i+angle_edge)])>>1;
            int mid_pix = *(bdata + (mid_y * iw) + mid_x);
            /* If the interior pixel value is not the same as the feature's... */
            if (mid_pix != feature_pix) {
                /* Reject the high-curvature minutia and return LQM_IGNORE. */
                return LQM_IGNORE;
            }

            /* Compute new direction based on line connecting adjusted feature */
            /* location and the midpoint in the feature's interior.            */
            int idir = LineToDirection(contour.xV[static_cast<std::size_t>(min_i)], contour.yV[static_cast<std::size_t>(min_i)], mid_x, mid_y, lqmParams.num_directions);

            /* Set minutia location to minimum theta position on the contour. */
            *oidir = idir;
            *ox_loc = contour.xV[static_cast<std::size_t>(min_i)];
            *oy_loc = contour.yV[static_cast<std::size_t>(min_i)];
            *ox_edge = contour.exV[static_cast<std::size_t>(min_i)];
            *oy_edge = contour.eyV[static_cast<std::size_t>(min_i)];

            /*Return normally. */
            return 0;
        }

        int GetHighCurvatureContour(Contour& ocontour,
                        const int half_contour,
                        const int x_loc, const int y_loc,
                        const int x_edge, const int y_edge,
                        unsigned char *bdata, const int iw, const int ih)
        {
            int ret;

            /* Compute maximum length of complete contour */
            /* (2 half contours + feature point).         */
            int max_contour = (half_contour<<1) + 1;

            /* Get 1st half contour with clockwise neighbor trace. */
            Contour half1;
            ret = TraceContour(half1,
                half_contour, x_loc, y_loc, x_loc, y_loc, x_edge, y_edge,
                LQM_SCAN_CLOCKWISE, bdata, iw, ih);
            if (ret) {

                /* If trace was not possible ... */
                if (ret == LQM_IGNORE) {
                    /* Return with contour length equal to 0. */
                    return 0;
                }

                /* If 1st half contour forms a loop ... */
                if (ret == LQM_LOOP_FOUND) {
                    /* Need to reverse the 1st half contour so that the points are    */
                    /* in consistent order.                                           */
                    /* We need to add the original feature point to the list, so      */
                    /* allocate new contour list with reserved length of one plus     */
                    /* length of 1st half contour                                     */
                    ocontour.Reserve(half1.ncontour + 1);

                    /* We have the new contour allocated, so store the original feature point. */
                    ocontour.AddPoint(x_loc, y_loc, x_edge, y_edge);

                    /* Now store the first half contour in reverse order. */
                    for (int j = half1.ncontour - 1; j >= 0; --j) {
                        ocontour.AddPoint(half1.xV[static_cast<std::size_t>(j)], half1.yV[static_cast<std::size_t>(j)], half1.exV[static_cast<std::size_t>(j)], half1.eyV[static_cast<std::size_t>(j)]);
                    }

                    /* Return LOOP_FOUND for further processing. */
                    return LQM_LOOP_FOUND;
                }

                /* Should be unreachable. Throw exception if reached here */
                throw std::logic_error(std::string("GetHighCurvatureContour: Received an unexpected return value from TraceContour: ") + std::to_string(ret));
            }

            /* If 1st half contour not complete ... */
            if (half1.ncontour < half_contour) {
                /* Return with contour length equal to 0. */
                return 0;
            }

            /* Otherwise, we have a complete 1st half contour...           */
            /* Get 2nd half contour with counter-clockwise neighbor trace. */
            /* Use the last point from the first contour trace as the      */
            /* point to test for a loop when tracing the second contour.   */
            Contour half2;
            ret = TraceContour(half2,
                half_contour, half1.xV[static_cast<std::size_t>(half1.ncontour - 1)], half1.yV[static_cast<std::size_t>(half1.ncontour - 1)],
                x_loc, y_loc, x_edge, y_edge,
                LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih);
            if (ret) {

                /* If 2nd trace was not possible ... */
                if (ret == LQM_IGNORE) {
                    /* Return with contour length equal to 0. */
                    return 0;
                }

                /* If non-zero return code is NOT LQM_LOOP_FOUND, then system error ... */
                if (ret != LQM_LOOP_FOUND) {
                    throw std::logic_error(std::string("GetHighCurvatureContour: Received an unexpected return value from TraceContour: ") + std::to_string(ret));
                }
            }

            /* If 2nd half NOT a loop AND not complete ... */
            if ((ret != LQM_LOOP_FOUND) && (half2.ncontour < half_contour)) {
                /* Return, with nothing allocated and contour length equal to 0. */
                return 0;
            }

            /* Otherwise we have a full 1st half contour and a 2nd half contour */
            /* that is either a loop or complete.  In either case we need to    */
            /* concatenate the two half contours into one longer contour.       */

            /* Allocate output contour list.  Go ahead and allocate the    */
            /* "max_contour" amount even though the resulting contour will */
            /* likely be shorter if it forms a loop.                       */
            ocontour.Reserve(max_contour);

            /* Copy 1st half contour into output contour buffers.      */
            /* This contour was collected clockwise, so it's points    */
            /* are entered in reverse order of the trace.  The result  */
            /* is the first point in the output contour if farthest    */
            /* from the starting feature point.                        */
            for (int j = half1.ncontour - 1; j >= 0; --j) {
                ocontour.AddPoint(half1.xV[static_cast<std::size_t>(j)], half1.yV[static_cast<std::size_t>(j)], half1.exV[static_cast<std::size_t>(j)], half1.eyV[static_cast<std::size_t>(j)]);
            }

            /* Next, store starting feature point into output contour buffers. */
            ocontour.AddPoint(x_loc, y_loc, x_edge, y_edge);

            /* Now, append 2nd half contour to permanent contour buffers.  */
            for (int j = 0; j < half2.ncontour; ++j) {
                ocontour.AddPoint(half2.xV[static_cast<std::size_t>(j)], half2.yV[static_cast<std::size_t>(j)], half2.exV[static_cast<std::size_t>(j)], half2.eyV[static_cast<std::size_t>(j)]);
            }

            /* Return the resulting return code form the 2nd call to TraceContour  */
            /* (the value will either be 0 or LQM_LOOP_FOUND).                     */
            // return ret; //todo: move this to second branch for fixes to bugs in original code. Reproduce original behavior for now
            return 0; // Original ret overwritten by allocate_contour()
        }

        int TraceContour(Contour& ocontour,
                        const int max_len, const int x_loop, const int y_loop,
                        const int x_loc, const int y_loc,
                        const int x_edge, const int y_edge,
                        const int scan_clock,
                        unsigned char *bdata, const int iw, const int ih)
        {
            /* Check to make sure that the feature and edge values are opposite. */
            if (*(bdata+(y_loc*iw)+x_loc) == *(bdata+(y_edge*iw)+x_edge)) {
                /* If not opposite, then the trace will not work, so return IGNORE. */
                return LQM_IGNORE;
            }

            ocontour.Reserve(max_len);

            /* Set up for finding first contour pixel. */
            int cur_x_loc = x_loc;
            int cur_y_loc = y_loc;
            int cur_x_edge = x_edge;
            int cur_y_edge = y_edge;

            /* Foreach pixel to be collected on the feature's contour... */
            for (int i = 0; i < max_len; ++i) {
                /* Find the next contour pixel. */
                int next_x_loc, next_y_loc, next_x_edge, next_y_edge;
                if (NextContourPixel(next_x_loc, next_y_loc,
                                        next_x_edge, next_y_edge,
                                        cur_x_loc, cur_y_loc,
                                        cur_x_edge, cur_y_edge,
                                        scan_clock, bdata, iw, ih)) {

                    /* If we trace back around to the specified starting */
                    /* feature location...                               */
                    if ((next_x_loc == x_loop) && (next_y_loc == y_loop)) {
                        /* Then we have found a loop, so return what we */
                        /* have traced to this point.                   */
                        return LQM_LOOP_FOUND;
                    }

                    /* Otherwise, we found another point on our feature's contour, */
                    /* so store the new contour point.                             */
                    ocontour.AddPoint(next_x_loc, next_y_loc, next_x_edge, next_y_edge);

                    /* Set up for finding next contour pixel. */
                    cur_x_loc = next_x_loc;
                    cur_y_loc = next_y_loc;
                    cur_x_edge = next_x_edge;
                    cur_y_edge = next_y_edge;
                }
                /* Otherwise, no new contour point found ... */
                else {
                    /* So, stop short and return normally with what we have   */
                    return 0;
                }
            }

            /* If we get here, we successfully found the maximum points we    */
            /* were looking for on the feature contour, so return normally.   */
            return 0;
        }

        int NextContourPixel(int& next_x_loc, int& next_y_loc,
                        int& next_x_edge, int& next_y_edge,
                        const int cur_x_loc, const int cur_y_loc,
                        const int cur_x_edge, const int cur_y_edge,
                        const int scan_clock,
                        unsigned char *bdata, const int iw, const int ih)
        {
            /* Get the feature's pixel value. */
            int feature_pix = *(bdata + (cur_y_loc * iw) + cur_x_loc);
            /* Get the feature's edge pixel value. */
            int edge_pix = *(bdata + (cur_y_edge * iw) + cur_x_edge);

            /* Get the nieghbor position of the feature's edge pixel in relationship */
            /* to the feature's actual position.                                     */
            /* REMEBER: The feature's position is always interior and on a ridge     */
            /* ending (black pixel) or (for bifurcations) on a valley ending (white  */
            /* pixel).  The feature's edge pixel is an adjacent pixel to the feature */
            /* pixel that is exterior to the ridge or valley ending and opposite in  */
            /* pixel value.                                                          */
            int nbr_i = StartScanNeighbor(cur_x_loc, cur_y_loc, cur_x_edge, cur_y_edge);

            /* Set current neighbor scan pixel to the feature's edge pixel. */
            int cur_nbr_x = cur_x_edge;
            int cur_nbr_y = cur_y_edge;
            int cur_nbr_pix = edge_pix;

            /* Foreach pixel neighboring the feature pixel ... */
            for (int i = 0; i < 8; ++i) {

                /* Set current neighbor scan pixel to previous scan pixel. */
                int prev_nbr_x = cur_nbr_x;
                int prev_nbr_y = cur_nbr_y;
                int prev_nbr_pix = cur_nbr_pix;

                /* Bump pixel neighbor index clockwise or counter-clockwise. */
                nbr_i = NextScanNeighbor(nbr_i, scan_clock);

                /* Set current scan pixel to the new neighbor.                   */
                /* REMEMBER: the neighbors are being scanned around the original */
                /* feature point.                                                */
                cur_nbr_x = cur_x_loc + OpenLQM::NEIGH8_DX[nbr_i];
                cur_nbr_y = cur_y_loc + OpenLQM::NEIGH8_DY[nbr_i];

                /* If new neighbor is not within image boundaries... */
                if ( (cur_nbr_x < 0) || (cur_nbr_x >= iw) || (cur_nbr_y < 0) || (cur_nbr_y >= ih) ) {
                    /* Return (LQM_FALSE==>Failure) if neighbor out of bounds. */
                    return LQM_FALSE;
                }

                /* Get the new neighbor's pixel value. */
                cur_nbr_pix = *(bdata + (cur_nbr_y * iw) + cur_nbr_x);

                /* If the new neighbor's pixel value is the same as the feature's   */
                /* pixel value AND the previous neighbor's pixel value is the same  */
                /* as the features's edge, then we have "likely" found our next     */
                /* contour pixel.                                                   */
                if ((cur_nbr_pix == feature_pix) && (prev_nbr_pix == edge_pix)) {

                    /* Check to see if current neighbor is on the corner of the */
                    /* neighborhood, and if so, test to see if it is "exposed". */
                    /* The neighborhood corners have odd neighbor indicies.     */
                    if (nbr_i % 2) {
                        /* To do this, look ahead one more neighbor pixel. */
                        int ni = NextScanNeighbor(nbr_i, scan_clock);
                        int nx = cur_x_loc + OpenLQM::NEIGH8_DX[ni];
                        int ny = cur_y_loc + OpenLQM::NEIGH8_DY[ni];
                        /* If new neighbor is not within image boundaries... */
                        if ( (nx < 0) || (nx >= iw) || (ny < 0) || (ny >= ih) ) {
                            /* Return (LQM_FALSE==>Failure) if neighbor out of bounds. */
                            return LQM_FALSE;
                        }
                        int npix = *(bdata + (ny * iw) + nx);

                        /* If the next neighbor's value is also the same as the */
                        /* feature's pixel, then corner is NOT exposed...       */
                        if (npix == feature_pix) {
                            /* Assign the current neighbor pair to the output pointers. */
                            next_x_loc = cur_nbr_x;
                            next_y_loc = cur_nbr_y;
                            next_x_edge = prev_nbr_x;
                            next_y_edge = prev_nbr_y;
                            /* Return LQM_TRUE==>Success. */
                            return LQM_TRUE;
                        }
                        /* Otherwise, corner pixel is "exposed" so skip it. */
                        else {
                            /* Skip current corner neighbor by resetting it to the      */
                            /* next neighbor, which upon the iteration will immediately */
                            /* become the previous neighbor.                            */
                            cur_nbr_x = nx;
                            cur_nbr_y = ny;
                            cur_nbr_pix = npix;
                            /* Advance neighbor index. */
                            nbr_i = ni;
                            /* Advance neighbor count. */
                            ++i;
                        }
                    }
                    /* Otherwise, current neighbor is not a corner ... */
                    else {
                        /* Assign the current neighbor pair to the output pointers. */
                        next_x_loc = cur_nbr_x;
                        next_y_loc = cur_nbr_y;
                        next_x_edge = prev_nbr_x;
                        next_y_edge = prev_nbr_y;
                        /* Return LQM_TRUE==>Success. */
                        return LQM_TRUE;
                    }
                }
            }

            /* If we get here, then we did not find the next contour pixel */
            /* within the 8 neighbors of the current feature pixel so      */
            /* return (FALSE==>Failure).                                   */
            /* NOTE: This must mean we found a single isolated pixel.      */
            /*       Perhaps this should be filled?                        */
            return LQM_FALSE;
        }

        int StartScanNeighbor(const int x_prev, const int y_prev, const int x_next, const int y_next)
        {
            if ((x_prev == x_next) && (y_next > y_prev)) {
                return LQM_SOUTH;
            }
            else if ((x_prev == x_next) && (y_next < y_prev)) {
                return LQM_NORTH;
            }
            else if ((x_next > x_prev) && (y_prev == y_next)) {
                return LQM_EAST;
            }
            else if ((x_next < x_prev) && (y_prev == y_next)) {
                return LQM_WEST;
            }

            throw std::invalid_argument("StartScanNeighbor: Invalid combination of next and prev points received");
        }

        int NextScanNeighbor(const int nbr_i, const int scan_clock)
        {
            /* If scanning neighbors clockwise, then advance one neighbor clockwise. */

            /* Otherwise, advance one neighbor counter-clockwise. */
            /* There are 8 pixels in the neighborhood, so to      */
            /* decrement with wrapping from 0 around to 7, add    */
            /* 7 to the neighbor index and mod with 8.            */

            return (nbr_i + (scan_clock == LQM_SCAN_CLOCKWISE ? 1 : 7)) % 8;
        }

        int IsLoopClockwise(const int *contour_x, const int *contour_y, const int ncontour, const int default_ret)
        {
            /* Derive chain code from contour points. */
            std::vector<int> chain;
            ChainCodeLoop(chain, contour_x, contour_y, ncontour);

            /* If chain is empty... */
            if (chain.size() == 0) {
                /* There wasn't enough contour points to tell, so return the   */
                /* the default return value.                                   */
                return default_ret;
            }

            /* If the chain code for contour is clockwise ... pass default return   */
            /* value on to this routine to correctly handle the case where we can't */
            /* tell the direction of the chain code.                                */
            return IsChainClockwise(chain, default_ret);
        }

        void ChainCodeLoop(std::vector<int>& ochain, const int *contour_x, const int *contour_y, const int ncontour)
        {
            /* If we don't have at least 3 points in the contour ... */
            if (ncontour <= 3) {
                /* Then we don't have a loop, so return */
                return;
            }

            /* Allocate chain code vector.  It will be the same length as the */
            /* number of points in the contour.  There will be one chain code */
            /* between each point on the contour including a code between the */
            /* last to the first point on the contour (completing the loop).  */
            ochain.reserve(static_cast<std::size_t>(ncontour));

            /* For each neighboring point in the list (with "i" pointing to the */
            /* previous neighbor and "j" pointing to the next neighbor...       */
            int i = 0, dx, dy;
            for (int j = 1; i < ncontour - 1; ++i, ++j) {
                /* Compute delta in X between neighbors. */
                dx = contour_x[j] - contour_x[i];
                /* Compute delta in Y between neighbors. */
                dy = contour_y[j] - contour_y[i];
                /* Derive chain code index from neighbor deltas.                  */
                /* The deltas are on the range [-1..1], so to use them as indices */
                /* into the code list, they must first be incremented by one.     */
                ochain.push_back(*(OpenLQM::NEIGH8_CHAINCODES + ((dy+1) * OpenLQM::NEIGH8_DIM) + dx + 1));
            }

            /* Now derive chain code between last and first points in the */
            /* contour list.                                              */
            dx = contour_x[0] - contour_x[i];
            dy = contour_y[0] - contour_y[i];
            ochain.push_back(*(OpenLQM::NEIGH8_CHAINCODES + ((dy + 1) * OpenLQM::NEIGH8_DIM) + dx + 1));
        }

        int IsChainClockwise(const std::vector<int>& chain, const int default_ret)
        {
            /* Initialize turn-accumulator to 0. */
            int sum = 0;

            /* Foreach neighboring code in chain, compute the difference in  */
            /* direction and accumulate.  Left-hand turns increment, whereas */
            /* right-hand decrement.                                         */
            int i = 0;
            for (int j = 1; i < (static_cast<int>(chain.size())) - 1; ++i, ++j) {
                /* Compute delta in neighbor direction. */
                int d = chain[static_cast<std::size_t>(j)] - chain[static_cast<std::size_t>(i)];
                /* Make the delta the "inner" distance. */
                /* If delta >= 4, for example if chain_i==2 and chain_j==7 (which   */
                /* means the contour went from a step up to step down-to-the-right) */
                /* then 5=(7-2) which is >=4, so -3=(5-8) which means that the      */
                /* change in direction is a righ-hand turn of 3 units).             */
                if (d >= 4) {
                    d -= 8;
                }
                /* If delta <= -4, for example if chain_i==7 and chain_j==2 (which  */
                /* means the contour went from a step down-to-the-right to step up) */
                /* then -5=(2-7) which is <=-4, so 3=(-5+8) which means that the    */
                /* change in direction is a left-hand turn of 3 units).             */
                else if (d <= -4) {
                    d += 8;
                }

                /* The delta direction is then accumulated. */
                sum += d;
            }

            /* Now we need to add in the final delta direction between the last */
            /* and first codes in the chain.                                    */
            int d = chain[0] - chain[static_cast<std::size_t>(i)];
            if (d >= 4) {
                d -= 8;
            }
            else if (d <= -4) {
                d += 8;
            }
            sum += d;

            /* If the final turn_accumulator == 0, then we CAN'T TELL the       */
            /* direction of the chain code, so return the default return value. */
            if (sum == 0) {
                return default_ret;
            }
            /* Otherwise, if the final turn-accumulator is positive ... */
            else if (sum > 0) {
                /* Then we had a greater amount of left-hand turns than right-hand     */
                /* turns, so the chain is in COUNTER-CLOCKWISE order, so return FALSE. */
                return LQM_FALSE;
            }
            /* Otherwise, the final turn-accumulator is negative ... */
            else {
                /* So we had a greater amount of right-hand turns than left-hand  */
                /* turns, so the chain is in CLOCKWISE order, so return LQM_TRUE.     */
                return LQM_TRUE;
            }
        }

        void ProcessMinutiaeLoop(Minutiae& minutiae,
                    const Contour& contour,
                    unsigned char *bdata, const int iw, const int ih,
                    int *plow_flow_map, const LQMParams& lqmParams)
        {
            /* If contour is empty, then just return. */
            if(contour.ncontour <= 0) {
                return;
            }

            /* If loop is large enough ... */
            if (contour.ncontour > lqmParams.min_loop_len) {
                /* Get pixel value of feature's interior. */
                int feature_pix = *(bdata + (contour.yV[0] * iw) + contour.xV[0]);

                /* Get the aspect dimensions of the loop in units of */
                /* squared distance.                                 */
                double min_dist, max_dist;
                int min_fr, max_fr, min_to, max_to;
                GetLoopAspect(&min_fr, &min_to, &min_dist,
                                &max_fr, &max_to, &max_dist,
                                contour);

                /* If loop passes aspect ratio tests ... loop is sufficiently  */
                /* narrow or elongated ...                                     */
                if ((min_dist < lqmParams.min_loop_aspect_dist) || ((max_dist/min_dist) >= lqmParams.min_loop_aspect_ratio)) {
                    /* Update minutiae list with opposite points of max distance */
                    /* on the loop.                                             */

                    /* First, check if interior point has proper pixel value. */
                    int mid_x = (contour.xV[static_cast<std::size_t>(max_fr)] + contour.xV[static_cast<std::size_t>(max_to)]) >> 1;
                    int mid_y = (contour.yV[static_cast<std::size_t>(max_fr)] + contour.yV[static_cast<std::size_t>(max_to)]) >> 1;
                    int mid_pix = *(bdata + (mid_y * iw) + mid_x);
                    /* If interior point is the same as the feature... */
                    if (mid_pix == feature_pix) {
                        /* 1. Treat maximum distance point as a potential minutia. */

                        /* Compute direction from maximum loop point to its */
                        /* opposite point.                                  */
                        int idir = LineToDirection(contour.xV[static_cast<std::size_t>(max_fr)], contour.yV[static_cast<std::size_t>(max_fr)],
                                            contour.xV[static_cast<std::size_t>(max_to)], contour.yV[static_cast<std::size_t>(max_to)],
                                            lqmParams.num_directions);
                        /* Get type of minutia: BIFURCATION or RIDGE_ENDING. */
                        int type = static_cast<bool>(feature_pix);
                        /* Determine if minutia is appearing or disappearing. */
                        int appearing = IsMinutiaAppearing(contour.xV[static_cast<std::size_t>(max_fr)], contour.yV[static_cast<std::size_t>(max_fr)], contour.exV[static_cast<std::size_t>(max_fr)], contour.eyV[static_cast<std::size_t>(max_fr)]);
                        /* Is the new point in a LOW RIDGE FLOW block? */
                        int fmapval = *(plow_flow_map+(contour.yV[static_cast<std::size_t>(max_fr)] * iw) + contour.xV[static_cast<std::size_t>(max_fr)]);

                        double reliability;
                        /* If current minutia is in a LOW RIDGE FLOW block ... */
                        if (fmapval) {
                            reliability = OpenLQM::Presets::MinutiaeReliability::MEDIUM;
                        }
                        else {
                            /* Otherwise, minutia is in a reliable block. */
                            reliability = OpenLQM::Presets::MinutiaeReliability::HIGH;
                        }

                        {
                            /* Create new minutia object. */
                            Minutia minutia(
                                contour.xV[static_cast<std::size_t>(max_fr)], contour.yV[static_cast<std::size_t>(max_fr)],
                                contour.exV[static_cast<std::size_t>(max_fr)], contour.eyV[static_cast<std::size_t>(max_fr)],
                                idir, reliability,
                                type, appearing, OpenLQM::Constants::Minutia::LOOP_ID
                            );

                            /* Update the minutiae list with potential new minutia.  */
                            UpdateMinutiae(minutiae, minutia, bdata, iw, ih, lqmParams);
                        }

                        /* 2. Treat point opposite of maximum distance point as */
                        /*    a potential minutia.                              */

                        /* Flip the direction 180 degrees. Make sure new direction */
                        /* is on the range [0..(ndirsX2)].                         */
                        idir += lqmParams.num_directions;
                        idir %= (lqmParams.num_directions << 1);

                        /* The type of minutia will stay the same. */

                        /* Determine if minutia is appearing or disappearing. */
                        appearing = IsMinutiaAppearing(contour.xV[static_cast<std::size_t>(max_to)], contour.yV[static_cast<std::size_t>(max_to)], contour.exV[static_cast<std::size_t>(max_to)], contour.eyV[static_cast<std::size_t>(max_to)]);

                        /* Is the new point in a LOW RIDGE FLOW block? */
                        fmapval = *(plow_flow_map+(contour.yV[static_cast<std::size_t>(max_to)]*iw)+
                                                    contour.xV[static_cast<std::size_t>(max_to)]);

                        /* If current minutia is in a LOW RIDGE FLOW block ... */
                        if (fmapval) {
                            reliability = OpenLQM::Presets::MinutiaeReliability::MEDIUM;
                        }
                        else {
                            /* Otherwise, minutia is in a reliable block. */
                            reliability = OpenLQM::Presets::MinutiaeReliability::HIGH;
                        }

                        {
                            /* Create new minutia object. */
                            Minutia minutia(
                                contour.xV[static_cast<std::size_t>(max_to)], contour.yV[static_cast<std::size_t>(max_to)],
                                contour.exV[static_cast<std::size_t>(max_to)], contour.eyV[static_cast<std::size_t>(max_to)],
                                idir, reliability,
                                type, appearing, OpenLQM::Constants::Minutia::LOOP_ID
                            );

                            /* Update the minutiae list with potential new minutia. */
                            /* NOTE: Deliberately using version one of this routine. */
                            UpdateMinutiae(minutiae, minutia, bdata, iw, ih, lqmParams);
                        }

                        /* Done successfully processing this loop, so return normally. */
                        return;
                    } /* Otherwise, loop interior has problems. */
                } /* Otherwise, loop is not the right shape for minutiae. */
            } /* Otherwise, loop's perimeter is too small for minutiae. */

            /* If we get here, we have a loop that is assumed to not contain */
            /* minutiae, so remove the loop from the image.                   */
            FillLoop(contour.xV.data(), contour.yV.data(), contour.ncontour, bdata, iw);
        }

        void GetLoopAspect(int *omin_fr, int *omin_to, double *omin_dist,
                    int *omax_fr, int *omax_to, double *omax_dist,
                    const Contour& contour)
        {
            /* Compute half the perimeter of the loop. */
            int halfway = contour.ncontour>>1;

            /* Take opposite points on the contour and walk half way    */
            /* around the loop.                                         */
            int i = 0;
            int j = halfway;
            /* Compute squared distance between opposite points on loop. */
            double dist = SquaredDistance(contour.xV[static_cast<std::size_t>(i)], contour.yV[static_cast<std::size_t>(i)], contour.xV[static_cast<std::size_t>(j)], contour.yV[static_cast<std::size_t>(j)]);

            /* Initialize running minimum and maximum distances along loop. */
            double min_dist = dist;
            int min_i = i;
            int min_j = j;
            double max_dist = dist;
            int max_i = i;
            int max_j = j;
            /* Bump to next pair of opposite points. */
            ++i;
            /* Make sure j wraps around end of list. */
            ++j;
            j %= contour.ncontour;

            /* If the loop is of even length, then we only need to walk half */
            /* way around as the other half will be exactly redundant.  If   */
            /* the loop is of odd length, then the second half will not be   */
            /* be exactly redundant and the difference "may" be meaningful.  */
            /* If execution speed is an issue, then probably get away with   */
            /* walking only the fist half of the loop under ALL conditions.  */

            int limit;
            /* If loop has odd length ... */
            if (contour.ncontour % 2) {
                /* Walk the loop's entire perimeter. */
                limit = contour.ncontour;
            }
            /* Otherwise the loop has even length ... */
            else {
                /* Only walk half the perimeter. */
                limit = halfway;
            }

            /* While we have not reached our perimeter limit ... */
            while (i < limit) {
                /* Compute squared distance between opposite points on loop. */
                dist = SquaredDistance(contour.xV[static_cast<std::size_t>(i)], contour.yV[static_cast<std::size_t>(i)], contour.xV[static_cast<std::size_t>(j)], contour.yV[static_cast<std::size_t>(j)]);
                /* Check the running minimum and maximum distances. */
                if(dist < min_dist){
                    min_dist = dist;
                    min_i = i;
                    min_j = j;
                }
                if(dist > max_dist){
                    max_dist = dist;
                    max_i = i;
                    max_j = j;
                }
                /* Bump to next pair of opposite points. */
                ++i;
                /* Make sure j wraps around end of list. */
                ++j;
                j %= contour.ncontour;
            }

            /* Assign minimum and maximum distances to output pointers. */
            *omin_fr = min_i;
            *omin_to = min_j;
            *omin_dist = min_dist;
            *omax_fr = max_i;
            *omax_to = max_j;
            *omax_dist = max_dist;
        }

        int LineToDirection(const int fx, const int fy, const int tx, const int ty, const int ndirs)
        {
            /* Compute angle to line connecting the 2 points.             */
            /* Coordinates are swapped and order of points reversed to    */
            /* account for 0 direction is vertical and positive direction */
            /* is clockwise.                                              */
            double theta = AngleToLine(ty, tx, fy, fx);

            /* Make sure the angle is positive. */
            theta += LQM_PI2;
            theta = fmod(theta, LQM_PI2);
            /* Convert from radians to integer direction on range [0..(ndirsX2)]. */
            /* Multiply radians by units/radian ((ndirsX2)/(2PI)), and you get    */
            /* angle in integer units.                                            */
            /* Compute number of directions on full circle. */
            int full_ndirs = ndirs<<1;
            /* Compute the radians to integer direction conversion factor. */
            double pi_factor = static_cast<double>(full_ndirs)/LQM_PI2;
            /* Convert radian angle to integer direction on full circle. */
            theta *= pi_factor;
            /* Need to truncate precision so that answers are consistent */
            /* on different computer architectures when rounding doubles. */
            theta = trunc_dbl_precision(theta);
            int idir = sround(theta);
            /* Make sure on range [0..(ndirsX2)]. */
            idir %= full_ndirs;

            /* Return the integer direction. */
            return idir;
        }

        double AngleToLine(const int fx, const int fy, const int tx, const int ty)
        {
            /* Compute slope of line connecting the 2 specified points. */
            double dy = static_cast<double>(fy - ty);
            double dx = static_cast<double>(tx - fx);

            double theta;
            /* If delta's are sufficiently small ... */
            if ((std::fabs(dx) < OpenLQM::Constants::Minutia::MIN_SLOPE_DELTA) && (std::fabs(dy) < OpenLQM::Constants::Minutia::MIN_SLOPE_DELTA)) {
                theta = 0.0;
            }
            /* Otherwise, compute angle to the line. */
            else {
                theta = atan2(dy, dx);
            }

            /* Return the compute angle in radians. */
            return theta;
        }

        int IsMinutiaAppearing(const int x_loc, const int y_loc, const int x_edge, const int y_edge)
        {
            /* Edge pixels will always be N,S,E,W of feature pixel. */

            /* 1. When scanning for feature's HORIZONTALLY... */
            /* If the edge is above the feature, then appearing. */
            if (x_edge < x_loc) {
                return LQM_APPEARING;
            }
            /* If the edge is below the feature, then disappearing. */
            if (x_edge > x_loc) {
                return LQM_DISAPPEARING;
            }

            /* 1. When scanning for feature's VERTICALLY... */
            /* If the edge is left of feature, then appearing. */
            if (y_edge < y_loc) {
                return LQM_APPEARING;
            }
            /* If the edge is right of feature, then disappearing. */
            if (y_edge > y_loc) {
                return LQM_DISAPPEARING;
            }

            /* Should never get here, but just in case. */
            throw std::invalid_argument("IsMinutiaAppearing: bad configuration of pixels");
        }

        int UpdateMinutiae(Minutiae& minutiae, Minutia& minutia,
                        unsigned char *bdata, const int iw, const int ih,
                        const LQMParams& lqmParams)
        {
            int i, dy, dx, delta_dir;

            /* Compute quarter of possible directions in a semi-circle */
            /* (ie. 45 degrees).                                       */
            int qtr_ndirs = lqmParams.num_directions >> 2;

            /* Compute number of directions in full circle. */
            int full_ndirs = lqmParams.num_directions<<1;

            /* Is the minutiae list empty? */
            if (minutiae.size() > 0) {
                /* Foreach minutia stored in the list... */
                for (i = 0; i < static_cast<int>(minutiae.size()); ++i) {
                    /* If x distance between new minutia and current list minutia */
                    /* are sufficiently close...                                 */
                    dx = std::abs(minutiae[static_cast<std::size_t>(i)].x - minutia.x);
                    if (dx < lqmParams.max_minutia_delta) {
                        /* If y distance between new minutia and current list minutia */
                        /* are sufficiently close...                                 */
                        dy = std::abs(minutiae[static_cast<std::size_t>(i)].y - minutia.y);
                        if (dy < lqmParams.max_minutia_delta) {
                            /* If new minutia and current list minutia are same type... */
                            if (minutiae[static_cast<std::size_t>(i)].type == minutia.type) {
                                /* Test to see if minutiae have similar directions. */
                                /* Take minimum of computed inner and outer        */
                                /* direction differences.                          */
                                delta_dir = std::abs(minutiae[static_cast<std::size_t>(i)].direction -
                                                minutia.direction);
                                delta_dir = std::min(delta_dir, full_ndirs - delta_dir);
                                /* If directional difference is <= 45 degrees... */
                                if (delta_dir <= qtr_ndirs) {
                                    /* If new minutia and current list minutia share */
                                    /* the same point... */
                                    if ((dx==0) && (dy==0)) {
                                        /* Then the minutiae match, so don't add the new one */
                                        /* to the list.                                     */
                                        return LQM_IGNORE;
                                    }
                                    /* Othewise, check if they share the same contour. */
                                    /* Start by searching "max_minutia_delta" steps    */
                                    /* clockwise.                                      */
                                    /* If new minutia point found on contour...        */
                                    if (SearchContour(minutia.x, minutia.y,
                                            lqmParams.max_minutia_delta,
                                            minutiae[static_cast<std::size_t>(i)].x, minutiae[static_cast<std::size_t>(i)].y,
                                            minutiae[static_cast<std::size_t>(i)].ex, minutiae[static_cast<std::size_t>(i)].ey,
                                            LQM_SCAN_CLOCKWISE, bdata, iw, ih)) {
                                        /* Consider the new minutia to be the same as the */
                                        /* current list minutia, so don't add the new one */
                                        /* to the list.                                   */
                                        return LQM_IGNORE;
                                    }
                                    /* Now search "max_minutia_delta" steps counter-  */
                                    /* clockwise along contour.                       */
                                    /* If new minutia point found on contour...       */
                                    if (SearchContour(minutia.x, minutia.y,
                                            lqmParams.max_minutia_delta,
                                            minutiae[static_cast<std::size_t>(i)].x, minutiae[static_cast<std::size_t>(i)].y,
                                            minutiae[static_cast<std::size_t>(i)].ex, minutiae[static_cast<std::size_t>(i)].ey,
                                            LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih)) {
                                        /* Consider the new minutia to be the same as the */
                                        /* current list minutia, so don't add the new one */
                                        /* to the list.                                   */
                                        return LQM_IGNORE;
                                    }

                                    /* Otherwise, new minutia and current list minutia do */
                                    /* not share the same contour, so although they are   */
                                    /* similar in type and location, treat them as 2      */
                                    /* different minutia.                                 */

                                } /* Otherwise, directions are too different. */
                            } /* Otherwise, minutiae are different type. */
                        } /* Otherwise, minutiae too far apart in Y. */
                    } /* Otherwise, minutiae too far apart in X. */
                } /* End FOR minutia in list. */
            } /* Otherwise, minutiae list is empty. */

            /* Otherwise, assume new minutia is not in the list, so add it. */
            minutiae.push_back(std::move(minutia));

            /* New minutia was successfully added to the list. */
            /* Return normally. */
            return 0;
        }

        int SearchContour(const int x_search, const int y_search,
                        const int search_len,
                        const int x_loc, const int y_loc,
                        const int x_edge, const int y_edge,
                        const int scan_clock,
                        unsigned char *bdata, const int iw, const int ih)
        {
            int next_x_loc, next_y_loc;
            int next_x_edge, next_y_edge;

            /* Set up for finding first contour pixel. */
            int cur_x_loc = x_loc;
            int cur_y_loc = y_loc;
            int cur_x_edge = x_edge;
            int cur_y_edge = y_edge;

            /* Foreach point to be collected on the feature's contour... */
            for (int i = 0; i < search_len; ++i) {
                /* Find the next contour pixel. */
                if (NextContourPixel(next_x_loc, next_y_loc,
                                        next_x_edge, next_y_edge,
                                        cur_x_loc, cur_y_loc,
                                        cur_x_edge, cur_y_edge,
                                        scan_clock, bdata, iw, ih)) {

                    /* If we find the point we are looking for on the contour... */
                    if ((next_x_loc == x_search) && (next_y_loc == y_search)) {
                        /* Then return FOUND. */
                        return LQM_FOUND;
                    }

                    /* Otherwise, set up for finding next contour pixel. */
                    cur_x_loc = next_x_loc;
                    cur_y_loc = next_y_loc;
                    cur_x_edge = next_x_edge;
                    cur_y_edge = next_y_edge;
                }
                /* Otherwise, no new contour point found ... */
                else {
                    /* So, stop searching, and return NOT_FOUND. */
                    return LQM_NOT_FOUND;
                }
            }

            /* If we get here, we successfully searched the maximum points */
            /* without finding our desired point, so return NOT_FOUND.     */
            return LQM_NOT_FOUND;
        }

        void FillLoop(const int *contour_x, const int *contour_y,
                        const int ncontour, unsigned char *bdata,
                        const int iw)
        {
            /* Create a shape structure from loop's contour. */
            Shape shape(contour_x, contour_y, ncontour);

            /* Get feature pixel value (the value on the interior of the loop */
            /* to be filled).                                                 */
            int feature_pix = *(bdata+(contour_y[0]*iw)+contour_x[0]);
            /* Now get edge pixel value (the value on the exterior of the loop    */
            /* to be used to filled the loop).  We can get this value by flipping */
            /* the feature pixel value.                                           */
            int edge_pix = (feature_pix ? 0 : 1);

            /* Foreach row in shape... */
            for (int i = 0; i < static_cast<int>(shape.rows.size()); ++i) {
                /* Get y-coord of current row in shape. */
                int y = shape.rows[static_cast<std::size_t>(i)].y;

                /* There should always be at least 1 contour points in the row.    */
                /* If there isn't, then something is wrong, so post a warning and  */
                /* just return.  This is mostly for debug purposes.                */
                if (shape.rows[static_cast<std::size_t>(i)].xs.size() < 1) {
                    // "WARNING : fill_loop : unexpected shape, preempting loop fill\n");
                    /* This is unexpected, but not fatal, so return normally. */
                    return;
                }

                /* Reset x index on row to the left-most contour point in the row. */
                int j = 0;
                /* Get first x-coord corresponding to the first contour point on row. */
                int x = shape.rows[static_cast<std::size_t>(i)].xs[static_cast<std::size_t>(j)];
                /* Fill the first contour point on the row. */
                *(bdata+(y*iw)+x) = static_cast<unsigned char>(edge_pix);
                /* Set the index of last contour point on row. */
                int lastj = (static_cast<int>(shape.rows[static_cast<std::size_t>(i)].xs.size())) - 1;
                /* While last contour point on row has not been processed... */
                while (j < lastj) {
                    /* On each interation, we have filled up to the current   */
                    /* contour point on the row pointed to by "j", and now we */
                    /* need to determine if we need to skip some edge pixels  */
                    /* caused by a concavity in the shape or not.             */

                    /* Get the next pixel value on the row just right of the     */
                    /* last contour point filled.  We know there are more points */
                    /* on the row because we haven't processed the last contour  */
                    /* point on the row yet.                                     */
                    ++x;
                    int next_pix = *(bdata+(y*iw)+x);

                    /* If the next pixel is the same value as loop's edge pixels ... */
                    if (next_pix == edge_pix) {
                        /* Then assume we have found a concavity and skip to next */
                        /* contour point on row.                                  */
                        j++;
                        /* Fill the new contour point because we know it is on the */
                        /* feature's contour.                                      */
                        x = shape.rows[static_cast<std::size_t>(i)].xs[static_cast<std::size_t>(j)];
                        *(bdata+(y*iw)+x) = static_cast<unsigned char>(edge_pix);

                        /* Now we are ready to loop again. */
                    }

                    /* Otherwise, fill from current pixel up through the next contour */
                    /* point to the right on the row.                                 */
                    else {
                        /* Bump to the next contour point to the right on row. */
                        ++j;
                        /* Set the destination x-coord to the next contour point   */
                        /* to the right on row.  Realize that this could be the    */
                        /* same pixel as the current x-coord if contour points are */
                        /* adjacent.                                               */
                        int nx = shape.rows[static_cast<std::size_t>(i)].xs[static_cast<std::size_t>(j)];

                        /* Fill between current x-coord and next contour point to the */
                        /* right on the row (including the new contour point).*/
                        FillPartialRow(edge_pix, x, nx, y, bdata, iw);
                    }

                    /* Once we are here we have filled the row up to (and including) */
                    /* the contour point currently pointed to by "j".                */
                    /* We are now ready to loop again.                               */

                } /* End WHILE */
            } /* End FOR */
        }

        void FillPartialRow(const int fill_pix, const int frx, const int tox,
                const int y, unsigned char *bdata, const int iw)
        {
            int x;
            unsigned char *bptr;

            /* Set pixel pointer to starting x-coord on current row. */
            bptr = bdata+(y*iw)+frx;

            /* Foreach pixel between starting and ending x-coord on row */
            /* (including the end points) ...                           */
            for (x = frx; x <= tox; ++x) {
                /* Set current pixel with fill pixel value. */
                *bptr = static_cast<unsigned char>(fill_pix);
                /* Bump to next pixel in the row. */
                ++bptr;
            }
        }

        int MinContourTheta(int *omin_i, double *omin_theta,
                            const int angle_edge,  const int *contour_x,
                            const int *contour_y, const int ncontour)
        {
            /* If the contour length is too short for processing... */
            if (ncontour < (angle_edge<<1)+1) {
                /* Return IGNORE. */
                return LQM_IGNORE;
            }

            /* Intialize running minimum values. */
            double min_theta = LQM_PI;
            /* Need to truncate precision so that answers are consistent   */
            /* on different computer architectures when comparing doubles. */
            min_theta = trunc_dbl_precision(min_theta);
            int min_i = -1;

            /* Set left angle point to first contour point. */
            int pleft = 0;
            /* Set center angle point to "angle_edge" points into contour. */
            int pcenter = angle_edge;
            /* Set right angle point to "angle_edge" points from pcenter. */
            int pright = pcenter + angle_edge;

            /* Loop until the right angle point exceeds the contour list. */
            while (pright < ncontour) {
                /* Compute angle to first edge line (connecting pcenter to pleft). */
                double theta1 = AngleToLine(contour_x[pcenter],contour_y[pcenter],
                                    contour_x[pleft],contour_y[pleft]);
                /* Compute angle to second edge line (connecting pcenter to pright). */
                double theta2 = AngleToLine(contour_x[pcenter],contour_y[pcenter],
                                    contour_x[pright],contour_y[pright]);

                /* Compute delta between angles accounting for an inner */
                /* and outer distance between the angles.               */
                double dtheta = std::fabs(theta2 - theta1);
                dtheta = std::min(dtheta, (LQM_PI*2.0) - dtheta);
                /* Need to truncate precision so that answers are consistent   */
                /* on different computer architectures when comparing doubles. */
                dtheta = trunc_dbl_precision(dtheta);

                /* Keep track of running minimum theta. */
                if(dtheta < min_theta){
                    min_i = pcenter;
                    min_theta = dtheta;
                }

                /* Bump to next points on contour. */
                pleft++;
                pcenter++;
                pright++;
            }

            /* If no minimum found (then contour is perfectly flat) so minimum */
            /* to center point on contour.                                     */
            if(min_i == -1){
                *omin_i = ncontour>>1;
                *omin_theta = min_theta;
            }
            else{
                /* Assign minimum theta information to output pointers. */
                *omin_i = min_i;
                *omin_theta = min_theta;
            }

            /* Return successfully. */
            return 0;
        }

        int GetLowCurvatureDirection(const int scan_dir, const int appearing, const int imapval, const int ndirs)
        {
            /* Start direction out with IMAP value. */
            int idir = imapval;

            /* NOTE!                                                           */
            /* The logic in this routine should hold whether for ridge endings */
            /* or for bifurcations.  The examples in the comments show ridge   */
            /* ending conditions only.                                         */

            /* CASE I : Ridge flow in Quadrant I; directions [0..8]             */
            if (imapval <= (ndirs>>1)) {
                /*    I.A: HORIZONTAL scan                                       */
                if (scan_dir == LQM_SCAN_HORIZONTAL) {
                    /*       I.A.1: Appearing Minutia                             */
                    if (appearing) {
                        /*        Ex.    0 0 0                                     */
                        /*               0 1 0                                     */
                        /*               ? ?                                       */
                        /*              Ridge flow is up and to the right, whereas */
                        /*              actual ridge is running down and to the    */
                        /*              left.                                      */
                        /*              Thus: HORIZONTAL : appearing   : should be */
                        /*                    OPPOSITE the ridge flow direction.   */
                        idir += ndirs;
                    }
                    // Otherwise:
                        /*       I.A.2: Disappearing Minutia                       */
                        /*           Ex.   ? ?                                     */
                        /*               0 1 0                                     */
                        /*               0 0 0                                     */
                        /*              Ridge flow is up and to the right, which   */
                        /*              should be SAME direction from which ridge  */
                        /*              is projecting.                             */
                        /*              Thus: HORIZONTAL : disappearing : should   */
                        /*                    be the same as ridge flow direction. */
                } /* End if HORIZONTAL scan */
                /* Otherwise:                                                    */
                /*    I.B: VERTICAL scan                                         */
                else{
                    /*       I.B.1: Disappearing Minutia                          */
                    if (!appearing) {
                        /*        Ex.    0 0                                       */
                        /*             ? 1 0                                       */
                        /*             ? 0 0                                       */
                        /*              Ridge flow is up and to the right, whereas */
                        /*              actual ridge is projecting down and to the */
                        /*              left.                                      */
                        /*              Thus: VERTICAL : disappearing : should be  */
                        /*                    OPPOSITE the ridge flow direction.   */
                        idir += ndirs;
                    }
                    // Otherwise:
                        /*       I.B.2: Appearing Minutia                          */
                        /*           Ex. 0 0 ?                                     */
                        /*               0 1 ?                                     */
                        /*               0 0                                       */
                        /*              Ridge flow is up and to the right, which   */
                        /*              should be SAME direction the ridge is      */
                        /*              running.                                   */
                        /*              Thus: VERTICAL : appearing   : should be   */
                        /*                    be the same as ridge flow direction. */
                } /* End else VERTICAL scan */
            } /* End if Quadrant I */

            // Otherwise:
            /* CASE II : Ridge flow in Quadrant II; directions [9..15]          */
            else {
                /*   II.A: HORIZONTAL scan                                       */
                if (scan_dir == LQM_SCAN_HORIZONTAL) {
                    /*      II.A.1: Disappearing Minutia                          */
                    if (!appearing) {
                        /*           Ex. ? ?                                       */
                        /*               0 1 0                                     */
                        /*               0 0 0                                     */
                        /*              Ridge flow is down and to the right,       */
                        /*              whereas actual ridge is running up and to  */
                        /*              the left.                                  */
                        /*              Thus: HORIZONTAL : disappearing : should   */
                        /*                    be OPPOSITE the ridge flow direction.*/
                        idir += ndirs;
                    }
                    // Otherwise:
                    /*      II.A.2: Appearing Minutia                             */
                        /*           Ex. 0 0 0                                     */
                        /*               0 1 0                                     */
                        /*                 ? ?                                     */
                        /*              Ridge flow is down and to the right, which */
                        /*              should be same direction from which ridge  */
                        /*              is projecting.                             */
                        /*              Thus: HORIZONTAL : appearing : should be   */
                        /*                    the SAME as ridge flow direction.    */
                } /* End if HORIZONTAL scan */
                /* Otherwise:                                                    */
                /*   II.B: VERTICAL scan                                         */
                else {
                    /*      II.B.1: Disappearing Minutia                          */
                    if (!appearing) {
                        /*           Ex. ? 0 0                                     */
                        /*               ? 1 0                                     */
                        /*                 0 0                                     */
                        /*              Ridge flow is down and to the right,       */
                        /*              whereas actual ridge is running up and to  */
                        /*              the left.                                  */
                        /*              Thus: VERTICAL : disappearing : should be  */
                        /*                    OPPOSITE the ridge flow direction.   */
                        idir += ndirs;
                    }
                    // Otherwise:
                        /*      II.B.2: Appearing Minutia                          */
                        /*           Ex. 0 0                                       */
                        /*               0 1 ?                                     */
                        /*               0 0 ?                                     */
                        /*              Ridge flow is down and to the right, which */
                        /*              should be same direction the ridge is      */
                        /*              projecting.                                */
                        /*              Thus: VERTICAL : appearing   : should be   */
                        /*                    be the SAME as ridge flow direction. */
                } /* End else VERTICAL scan */
            } /* End else Quadrant II */

            /* Return resulting direction on range [0..31]. */
            return idir;
        }

        int ScanAndUpdateMinutiae(Minutiae& minutiae, Minutia& minutia,
                        const int scan_dir, const int dmapval,
                        unsigned char *bdata, const int iw, const int ih,
                        const LQMParams& lqmParams)
        {
            /* Compute quarter of possible directions in a semi-circle */
            /* (ie. 45 degrees).                                       */
            int qtr_ndirs = lqmParams.num_directions >> 2;

            /* Compute number of directions in full circle. */
            int full_ndirs = lqmParams.num_directions << 1;

            /* Is the minutiae list empty? */
            if (minutiae.size() > 0) {
                /* Foreach minutia stored in the list (in reverse order) ... */
                for (int i = (static_cast<int>(minutiae.size())) - 1; i >= 0; --i) {
                    /* If x distance between new minutia and current list minutia */
                    /* are sufficiently close...                                 */
                    int dx = std::abs(minutiae[static_cast<std::size_t>(i)].x - minutia.x);
                    if (dx < lqmParams.max_minutia_delta) {
                        /* If y distance between new minutia and current list minutia */
                        /* are sufficiently close...                                 */
                        int dy = abs(minutiae[static_cast<std::size_t>(i)].y - minutia.y);
                        if (dy < lqmParams.max_minutia_delta) {
                            /* If new minutia and current list minutia are same type... */
                            if (minutiae[static_cast<std::size_t>(i)].type == minutia.type) {
                                /* Test to see if minutiae have similar directions. */
                                /* Take minimum of computed inner and outer        */
                                /* direction differences.                          */
                                int delta_dir = std::abs(minutiae[static_cast<std::size_t>(i)].direction - minutia.direction);
                                delta_dir = std::min(delta_dir, full_ndirs-delta_dir);
                                /* If directional difference is <= 45 degrees... */
                                if (delta_dir <= qtr_ndirs) {
                                    /* If new minutia and current list minutia share */
                                    /* the same point... */
                                    if ((dx==0) && (dy==0)) {
                                        /* Then the minutiae match, so don't add the new one */
                                        /* to the list.                                     */
                                        return(LQM_IGNORE);
                                    }
                                    /* Othewise, check if they share the same contour. */
                                    /* Start by searching "max_minutia_delta" steps    */
                                    /* clockwise.                                      */
                                    /* If new minutia point found on contour...        */
                                    if (SearchContour(minutia.x, minutia.y,
                                            lqmParams.max_minutia_delta,
                                            minutiae[static_cast<std::size_t>(i)].x, minutiae[static_cast<std::size_t>(i)].y,
                                            minutiae[static_cast<std::size_t>(i)].ex, minutiae[static_cast<std::size_t>(i)].ey,
                                            LQM_SCAN_CLOCKWISE, bdata, iw, ih) ||
                                        SearchContour(minutia.x, minutia.y,
                                            lqmParams.max_minutia_delta,
                                            minutiae[static_cast<std::size_t>(i)].x, minutiae[static_cast<std::size_t>(i)].y,
                                            minutiae[static_cast<std::size_t>(i)].ex, minutiae[static_cast<std::size_t>(i)].ey,
                                            LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih)) {
                                        /* If new minutia has VALID block direction ... */
                                        if (dmapval >= 0) {
                                            /* Derive feature scan direction compatible */
                                            /* with VALID direction.                    */
                                            int map_scan_dir = ChooseScanDirection(dmapval,
                                                                        lqmParams.num_directions);
                                            /* If map scan direction compatible with scan   */
                                            /* direction in which new minutia was found ... */
                                            if (map_scan_dir == scan_dir) {
                                                /* Then choose the new minutia over the one */
                                                /* currently in the list.                   */
                                                minutiae.erase(minutiae.begin() + i);
                                                /* Continue on ... */
                                            }
                                            else {
                                                /* Othersize, scan directions not compatible...*/
                                                /* so choose to keep the current minutia in    */
                                                /* the list and ignore the new one.            */
                                                return LQM_IGNORE;
                                            }
                                        }
                                        else {
                                            /* Otherwise, no reason to believe new minutia    */
                                            /* is any better than the current one in the list,*/
                                            /* so consider the new minutia to be the same as  */
                                            /* the current list minutia, and don't add the new*/
                                            /*  one to the list.                              */
                                            return LQM_IGNORE;
                                        }
                                    }

                                    /* Otherwise, new minutia and current list minutia do */
                                    /* not share the same contour, so although they are   */
                                    /* similar in type and location, treat them as 2      */
                                    /* different minutia.                                 */

                                } /* Otherwise, directions are too different. */
                            } /* Otherwise, minutiae are different type. */
                        } /* Otherwise, minutiae too far apart in Y. */
                    } /* Otherwise, minutiae too far apart in X. */
                } /* End FOR minutia in list. */
            } /* Otherwise, minutiae list is empty. */

            /* Otherwise, assume new minutia is not in the list, or those that */
            /* were close neighbors were selectively removed, so add it.       */
            minutiae.push_back(std::move(minutia));

            /* New minutia was successfully added to the list. */
            /* Return normally. */
            return 0;
        }

        int ChooseScanDirection(const int imapval, const int ndirs)
        {
            /* Compute quarter of directions in semi-circle. */
            int qtr_ndirs = ndirs>>2;

            /* If ridge flow in block is relatively vertical, then we want */
            /* to scan for minutia features in the opposite direction      */
            /* (ie. HORIZONTALLY).                                         */
            if ((imapval <= qtr_ndirs) || (imapval > (qtr_ndirs*3))) {
                return LQM_SCAN_HORIZONTAL;
            }
            /* Otherwise, ridge flow is realtively horizontal, and we want */
            /* to scan for minutia features in the opposite direction      */
            /* (ie. VERTICALLY).                                           */
            else {
                return LQM_SCAN_VERTICAL;
            }
        }

        void ScanForMinutiaeVertically(Minutiae& minutiae,
                        unsigned char *bdata, const int iw, const int ih,
                        int *pdirection_map, int *plow_flow_map, int *phigh_curve_map,
                        const LQMParams& lqmParams)
        {
            int possible[LQM_NFEATURES], nposs;

            /* Set scan region to entire image. */
            int sx = 0;
            int ex = iw;
            int sy = 0;
            int ey = ih;

            /* Start at first column in region. */
            int cx = sx;
            /* While second scan column not outside the right of the region ... */
            while (cx+1 < ex) {
                /* Start at beginning of new scan column in region. */
                int cy = sy;
                /* While not at end of region's current scan column. */
                while (cy < ey) {
                    /* Get pixel pair from current y position in current and next */
                    /* scan columns. */
                    unsigned char* p1ptr = bdata+(cy*iw)+cx;
                    unsigned char* p2ptr = p1ptr+1;
                    /* If scan pixel pair matches first pixel pair of */
                    /* 1 or more features... */
                    if (MatchFirstPair(*p1ptr, *p2ptr, possible, &nposs)) {
                        /* Bump forward to next scan pixel pair. */
                        cy++;
                        p1ptr+=iw;
                        p2ptr+=iw;
                        /* If not at end of region's current scan column... */
                        if (cy < ey) {
                            /* If scan pixel pair matches second pixel pair of */
                            /* 1 or more features... */
                            if (MatchSecondPair(*p1ptr, *p2ptr, possible, &nposs)) {
                                /* Store current y location. */
                                int y2 = cy;
                                /* Skip repeated pixel pairs. */
                                SkipRepeatedVerticalPair(&cy, ey, &p1ptr, &p2ptr,
                                                                iw);
                                /* If not at end of region's current scan column... */
                                if (cy < ey) {
                                    /* If scan pixel pair matches third pixel pair of */
                                    /* a single feature... */
                                    if (MatchThirdPair(*p1ptr, *p2ptr, possible, &nposs)) {
                                        /* Process detected minutia point. */
                                        ProcessVerticalScanMinutia(minutiae,
                                                        cx, cy, y2, possible[0],
                                                        bdata, iw, ih, pdirection_map,
                                                        plow_flow_map, phigh_curve_map,
                                                        lqmParams);
                                    }

                                    /* Set up to resume scan. */
                                    /* Test to see if 3rd pair can slide into 2nd pair. */
                                    /* The values of the 2nd pair MUST be different.    */
                                    /* If 3rd pair values are different ... */
                                    if (*p1ptr != *p2ptr) {
                                        /* Set next first pair to last of repeated */
                                        /* 2nd pairs, ie. back up one pair.        */
                                        --cy;
                                    }

                                    /* Otherwise, 3rd pair can't be a 2nd pair, so  */
                                    /* keep pointing to 3rd pair so that it is used */
                                    /* in the next first pair test.                 */

                                } /* Else, at end of current scan row. */
                            }

                            /* Otherwise, 2nd pair failed, so keep pointing to it */
                            /* so that it is used in the next first pair test.    */
                        } /* Else, at end of current scan column. */
                    }
                    /* Otherwise, 1st pair failed... */
                    else {
                        /* Bump forward to next pixel pair. */
                        ++cy;
                    }
                } /* While not at end of current scan column. */
                /* Bump forward to next scan column. */
                ++cx;
            } /* While not out of scan columns. */
        }

        void RemoveFalseMinutia(Minutiae& minutiae,
                unsigned char *bdata, const int iw, const int ih,
                int *direction_map, int *low_flow_map, int *high_curve_map,
                const int mw, const int mh, const LQMParams& lqmParams)
        {
            /* 1. Sort minutiae points top-to-bottom and left-to-right. */
            SortMinutiaeYX(minutiae, iw);

            /* 2. Remove minutiae on lakes (filled with white pixels) and        */
            /*    islands (filled with black pixels), both  defined by a pair of */
            /*    minutia points.                                                */
            RemoveIslandsAndLakes(minutiae, bdata, iw, ih, lqmParams);

            /* 3. Remove minutiae on holes in the binary image defined by a */
            /*    single point.                                             */
            RemoveHoles(minutiae, bdata, iw, ih, lqmParams);

            /* 4. Remove minutiae that point sufficiently close to a block with */
            /*    LQM_INVALID direction.                                            */
            RemovePointingInvalidBlock(minutiae, direction_map, mw, lqmParams);

            /* 5. Remove minutiae that are sufficiently close to a block with */
            /*    INVALID direction.                                          */
            RemoveNearInvalidBlocks(minutiae, direction_map, mw, mh, lqmParams);

            /* 6. Remove or adjust minutiae that reside on the side of a ridge */
            /*    or valley.                                                   */
            RemoveOrAdjustSideMinutiae(minutiae, bdata, iw, ih, direction_map, mw, lqmParams);

            /* 7. Remove minutiae that form a hook on the side of a ridge or valley. */
            RemoveHooks(minutiae, bdata, iw, ih, lqmParams);

            /* 8. Remove minutiae that are on opposite sides of an overlap. */
            RemoveOverlaps(minutiae, bdata, iw, lqmParams);

            /* 9. Remove minutiae that are "irregularly" shaped. */
            RemoveMalformations(minutiae, bdata, iw, ih, low_flow_map, mw, lqmParams);

            /* 10. Remove minutiae that form long, narrow, loops in the */
            /*     "unreliable" regions in the binary image.            */
            RemovePores(minutiae,  bdata, iw, ih,
                                        direction_map, low_flow_map, high_curve_map,
                                        mw, lqmParams);
        }

        void SortMinutiaeYX(Minutiae& minutiae, const int iw)
        {
            auto rank = [&](const Minutia& m) -> int {
                return m.y * iw + m.x;
            };

            std::sort(minutiae.begin(), minutiae.end(), [&](const Minutia& a, const Minutia& b) {
                return rank(a) < rank(b);
            });
        }

        void SortMinutiaeXY(Minutiae& minutiae, const int iw) {
            auto rank = [&](const Minutia& m) -> int {
                return m.x * iw + m.y;
            };

            std::sort(minutiae.begin(), minutiae.end(), [&](const Minutia& a, const Minutia& b) {
                return rank(a) < rank(b);
            });
        }

        void RemoveIslandsAndLakes(Minutiae& minutiae,
                            unsigned char *bdata, const int iw, const int ih,
                            const LQMParams& lqmParams)
        {
            int s, ret;

            int dist_thresh = lqmParams.max_rmtest_dist;
            int half_loop = lqmParams.max_half_loop;

            /* Allocate list of minutia indices that upon completion of testing  */
            /* should be removed from the minutiae lists. Initialized to false   */
            std::vector<bool> to_remove(minutiae.size(), false);

            /* Compute number directions in full circle. */
            int full_ndirs = lqmParams.num_directions << 1;
            /* Compute number of directions in 45=(180/4) degrees. */
            int qtr_ndirs = lqmParams.num_directions >> 2;

            /* Minimum allowable deltadir to consider joining minutia.               */
            /* (The closer the deltadir is to 180 degrees, the more likely the join. */
            /* When ndirs==16, then this value is 11=(3*4)-1 == 123.75 degrees.      */
            /* I chose to parameterize this threshold based on a fixed fraction of   */
            /* 'ndirs' rather than on passing in a parameter in degrees and doing    */
            /* the conversion.  I doubt the difference matters.                      */
            int min_deltadir = (3 * qtr_ndirs) - 1;

            /* Foreach primary (first) minutia (except for last one in list) ... */
            int f = 0;
            while (f < (static_cast<int>(minutiae.size())) - 1) {
                /* If current first minutia not previously set to be removed. */
                if (!to_remove[static_cast<std::size_t>(f)]) {
                    /* Set first minutia to temporary pointer. */
                    Minutia* minutia1 = &minutiae[static_cast<std::size_t>(f)];

                    /* Foreach secondary minutia to right of first minutia ... */
                    s = f+1;
                    while (s < static_cast<int>(minutiae.size())) {
                        /* Set second minutia to temporary pointer. */
                        Minutia* minutia2 = &minutiae[static_cast<std::size_t>(s)];

                        /* If the secondary minutia is desired type ... */
                        if (minutia2->type == minutia1->type) {
                            /* The binary image is potentially being edited during   */
                            /* each iteration of the secondary minutia loop,         */
                            /* therefore minutia pixel values may be changed.  We    */
                            /* need to catch these events by using the next 2 tests. */

                            /* If the first minutia's pixel has been previously */
                            /* changed...                                       */
                            if (*(bdata+(minutia1->y*iw)+minutia1->x) != minutia1->type) {
                                /* Then break out of secondary loop and skip to next */
                                /* first.                                            */
                                break;
                            }

                            /* If the second minutia's pixel has been previously */
                            /* changed...                                        */
                            if (*(bdata+(minutia2->y*iw)+minutia2->x) != minutia2->type) {
                                /* Set to remove second minutia. */
                                to_remove[static_cast<std::size_t>(s)] = true;
                            }

                            /* If the second minutia not previously set to be removed. */
                            if (!to_remove[static_cast<std::size_t>(s)]) {

                                /* Compute delta y between 1st & 2nd minutiae and test. */
                                int delta_y = minutia2->y - minutia1->y;
                                /* If delta y small enough (ex. <16 pixels)... */
                                if (delta_y <= dist_thresh) {
                                    /* Compute Euclidean distance between 1st & 2nd */
                                    /* mintuae.                                     */
                                    double dist = std::hypot(minutia1->x - minutia2->x, minutia1->y - minutia2->y);

                                    /* If distance is NOT too large (ex. <16 pixels)... */
                                    if (dist <= dist_thresh) {
                                        /* Compute "inner" difference between directions */
                                        /* on a full circle and test.                    */
                                        int inner_deltadir = OpenLQM::Core::ClosestDirDist(minutia1->direction, minutia2->direction, full_ndirs);
                                        if (inner_deltadir == LQM_INVALID_DIR) {
                                            throw std::invalid_argument("RemoveIslandsAndLakes: Invalid direction");
                                        }
                                        /* If the difference between dirs is large      */
                                        /* enough ...                                   */
                                        /* (the more 1st & 2nd point away from each     */
                                        /* other the more likely they should be joined) */
                                        if (inner_deltadir > min_deltadir) {
                                            /* Pair is the same type, so test to see */
                                            /* if both are on an island or lake.     */

                                            /* Check to see if pair on a loop of specified */
                                            /* half length (ex. 30 pixels) ...             */
                                            Contour loop;
                                            ret = OnIslandLake(loop,
                                                            *minutia1, *minutia2,
                                                            half_loop, bdata, iw, ih);
                                            /* If pair is on island/lake ... */
                                            if (ret == LQM_LOOP_FOUND) {
                                                /* Fill the loop. */
                                                FillLoop(loop.xV.data(), loop.yV.data(), loop.ncontour, bdata, iw);
                                                /* Set to remove first minutia. */
                                                to_remove[static_cast<std::size_t>(f)] = true;
                                                /* Set to remove second minutia. */
                                                to_remove[static_cast<std::size_t>(s)] = true;
                                            }
                                            /* If island/lake test IGNORED ... */
                                            else if (ret == LQM_IGNORE){
                                                /* Set to remove first minutia. */
                                                to_remove[static_cast<std::size_t>(f)] = true;
                                                /* Skip to next 1st minutia by breaking out */
                                                /* of inner secondary loop.                 */
                                                break;
                                            }
                                        }/* End deltadir test. */
                                    }/* End distance test. */
                                }
                                /* Otherwise, current 2nd too far below 1st, so skip to */
                                /* next 1st minutia.                                    */
                                else {
                                    /* Break out of inner secondary loop. */
                                    break;
                                }/* End delta-y test. */
                            }/* End if !to_remove[s] */
                        }/* End if 2nd not desired type */

                        /* Bump to next second minutia in minutiae list. */
                        s++;
                    }/* End secondary minutiae loop. */
                }/* Otherwise, first minutia already flagged to be removed. */

                /* Bump to next first minutia in minutiae list. */
                f++;
            }/* End primary minutiae loop. */

            /* Now remove all minutiae in list that have been flagged for removal. */
            /* NOTE: Need to remove the minutia from their lists in reverse       */
            /*       order, otherwise, indices will be off.                       */

            for (int i = (static_cast<int>(minutiae.size())) - 1; i >= 0; --i) {
                /* If the current minutia index is flagged for removal ... */
                if (to_remove[static_cast<std::size_t>(i)]) {
                    /* Remove the minutia from the minutiae list. */
                    minutiae.erase(minutiae.begin() + i);
                }
            }
        }

        int OnIslandLake(Contour& ocontour,
                        const Minutia& minutia1, const Minutia& minutia2,
                        const int max_half_loop,
                        unsigned char *bdata, const int iw, const int ih)
        {
            Contour contour1, contour2;

            /* Trace the contour of the feature starting at the 1st minutia point  */
            /* and stepping along up to the specified maximum number of steps or   */
            /* until 2nd mintuia point is encountered.                             */
            int ret = TraceContour(contour1, max_half_loop,
                                minutia2.x, minutia2.y, minutia1.x, minutia1.y,
                                minutia1.ex, minutia1.ey,
                                LQM_SCAN_CLOCKWISE, bdata, iw, ih);

            /* If trace was not possible, return IGNORE. */
            if (ret == LQM_IGNORE) {
                return ret;
            }

            /* If the trace encounters 2nd minutia point ... */
            if (ret == LQM_LOOP_FOUND) {

                /* Now, trace the contour of the feature starting at the 2nd minutia, */
                /* continuing to search for edge neighbors clockwise, and stepping    */
                /* along up to the specified maximum number of steps or until 1st     */
                /*  mintuia point is encountered.                                     */
                ret = TraceContour(contour2, max_half_loop,
                                    minutia1.x, minutia1.y, minutia2.x, minutia2.y,
                                    minutia2.ex, minutia2.ey,
                                    LQM_SCAN_CLOCKWISE, bdata, iw, ih);

                /* If trace was not possible, return IGNORE. */
                if (ret == LQM_IGNORE) {
                    return ret;
                }

                /* If the 2nd trace encounters 1st minutia point ... */
                if (ret == LQM_LOOP_FOUND) {
                    /* Combine the 2 half loop contours into one full loop. */

                    /* Compute loop length (including the minutia pair). */
                    int nloop = contour1.ncontour + contour2.ncontour + 2;

                    /* Allocate loop contour. */
                    ocontour.Reserve(nloop);

                    /* Store 1st minutia. */
                    ocontour.AddPoint(minutia1.x, minutia1.y, minutia1.ex, minutia1.ey);

                    /* Store first contour. */
                    for (int i = 0; i < contour1.ncontour; ++i) {
                        ocontour.AddPoint(contour1.xV[static_cast<std::size_t>(i)], contour1.yV[static_cast<std::size_t>(i)], contour1.exV[static_cast<std::size_t>(i)], contour1.eyV[static_cast<std::size_t>(i)]);
                    }
                    /* Store 2nd minutia. */
                    ocontour.AddPoint(minutia2.x, minutia2.y, minutia2.ex, minutia2.ey);
                    /* Store 2nd contour. */
                    for (int i = 0; i < contour2.ncontour; ++i) {
                        ocontour.AddPoint(contour2.xV[static_cast<std::size_t>(i)], contour2.yV[static_cast<std::size_t>(i)], contour2.exV[static_cast<std::size_t>(i)], contour2.eyV[static_cast<std::size_t>(i)]);
                    }

                    /* Then return that an island/lake WAS found (LQM_LOOP_FOUND). */
                    return LQM_LOOP_FOUND;
                }

                /* If the trace successfully followed 2nd minutia's contour, but   */
                /* did not encounter 1st minutia point within the specified number */
                /* of steps ...                                                    */
                if (ret == 0) {
                    /* Then return that an island/lake was NOT found (FALSE). */
                    return LQM_FALSE;
                }

                /* Otherwise, the 2nd trace had an error in following the contour ... */
                return ret;
            }

            /* If the 1st trace successfully followed 1st minutia's contour, but   */
            /* did not encounter the 2nd minutia point within the specified number */
            /* of steps ...                                                        */
            if (ret == 0) {
                /* Then return that an island/lake was NOT found (LQM_FALSE). */
                return LQM_FALSE;
            }

            /* Otherwise, the 1st trace had an error in following the contour ... */
            return ret;
        }

        void RemoveHoles(Minutiae& minutiae,
                        unsigned char *bdata, const int iw, const int ih,
                        const LQMParams& lqmParams)
        {
            int i = 0;
            /* Foreach minutia remaining in list ... */
            while (i < static_cast<int>(minutiae.size())) {
                /* Assign a temporary pointer. */
                Minutia& minutia = minutiae[static_cast<std::size_t>(i)];
                /* If current minutia is a bifurcation ... */
                if (minutia.type == LQM_BIFURCATION) {
                    /* Check to see if it is on a loop of specified length (ex. 15). */
                    int ret = MinutiaOnLoop(minutia, lqmParams.small_loop_len, bdata, iw, ih);
                    /* If minutia is on a loop ... or loop test IGNORED */
                    if ((ret == LQM_LOOP_FOUND) || (ret == LQM_IGNORE)) {
                        /* Then remove the minutia from list. */
                        minutiae.erase(minutiae.begin() + i);
                        /* No need to advance because next minutia has "slid" */
                        /* into position pointed to by 'i'.                   */
                    }
                    /* If the minutia is NOT on a loop... */
                    else if (ret == LQM_FALSE){
                        /* Simply advance to next minutia in the list. */
                        ++i;
                    }
                    /* Otherwise, an unknown value was returned! */
                    else {
                        throw std::logic_error(std::string("MinutiaOnLoop returned an unrecognized value: ") + std::to_string(ret));
                    }
                }
                /* Otherwise, the current minutia is a ridge-ending... */
                else {
                    /* Advance to next minutia in the list. */
                    ++i;
                }
            }
        }

        int MinutiaOnLoop(const Minutia& minutia, const int max_loop_len,
                    unsigned char *bdata, const int iw, const int ih)
        {
            int ret;
            Contour contour;

            /* Trace the contour of the feature starting at the minutia point  */
            /* and stepping along up to the specified maximum number of steps. */
            ret = TraceContour(contour, max_loop_len,
                                minutia.x, minutia.y, minutia.x, minutia.y,
                                minutia.ex, minutia.ey,
                                LQM_SCAN_CLOCKWISE, bdata, iw, ih);

            /* If trace was not possible ... */
            if (ret == LQM_IGNORE) {
                return ret;
            }

            /* If the trace completed a loop ... */
            if (ret == LQM_LOOP_FOUND) {
                return LQM_LOOP_FOUND;
            }

            /* If the trace successfully followed the minutia's contour, but did */
            /* not complete a loop within the specified number of steps ...      */
            if (ret == 0) {
                return LQM_FALSE;
            }

            /* Otherwise, TraceContour returned an unknown value */
            throw std::logic_error(std::string("TraceContour returned an unrecognized value: ") + std::to_string(ret));
        }

        void RemovePointingInvalidBlock(Minutiae& minutiae,
                                    int *direction_map, const int mw,
                                    const LQMParams& lqmParams)
        {
            /* Compute factor for converting integer directions to radians. */
            double pi_factor = LQM_PI / static_cast<double>(lqmParams.num_directions);

            int i = 0;
            /* Foreach minutia remaining in list ... */
            while (i < static_cast<int>(minutiae.size())) {
                /* Set temporary minutia pointer. */
                Minutia& minutia = minutiae[static_cast<std::size_t>(i)];
                /* Convert minutia's direction to radians. */
                double theta = minutia.direction * pi_factor;
                /* Compute translation offsets (ex. 6 pixels). */
                double dx = std::sin(theta) * static_cast<double>(lqmParams.trans_dir_pix);
                double dy = std::cos(theta) * static_cast<double>(lqmParams.trans_dir_pix);
                /* Need to truncate precision so that answers are consistent */
                /* on different computer architectures when rounding doubles. */
                dx = trunc_dbl_precision(dx);
                dy = trunc_dbl_precision(dy);
                int delta_x = sround(dx);
                int delta_y = sround(dy);
                /* Translate the minutia's coords. */
                int nx = minutia.x - delta_x;
                int ny = minutia.y + delta_y;
                /* Convert pixel coords to block coords. */
                int bx = static_cast<int>(nx / lqmParams.blocksize);
                int by = static_cast<int>(ny / lqmParams.blocksize);
                /* The translation could move the point out of image boundaries,    */
                /* and therefore the corresponding block coords can be out of       */
                /* map boundaries, so limit the block coords to within boundaries.  */

                /* Get corresponding block's ridge flow direction. */
                int dmapval = *(direction_map+(by*mw)+bx);

                /* If the block's direction is INVALID ... */
                if (dmapval == LQM_INVALID_DIR) {
                    /* Remove the minutia from the minutiae list. */
                    minutiae.erase(minutiae.begin() + i);
                    /* No need to advance because next minutia has slid into slot. */
                }
                else {
                    /* Advance to next minutia in list. */
                    ++i;
                }
            }
        }

        void RemoveNearInvalidBlocks(Minutiae& minutiae, int *direction_map, const int mw, const int mh, const LQMParams& lqmParams)
        {
            /* The next 2 lookup tables are indexed by 'ix' and 'iy'. */
            /* When a feature pixel lies within a 6-pixel margin of a */
            /* block, this routine examines neighboring blocks to     */
            /* determine appropriate actions.                         */
            /*    'ix' may take on values:                            */
            /*         0 == x-pixel coord in leftmost margin          */
            /*         1 == x-pixel coord in middle of block          */
            /*         2 == x-pixel coord in rightmost margin         */
            /*    'iy' may take on values:                            */
            /*         0 == y-pixel coord in topmost margin           */
            /*         1 == y-pixel coord in middle of block          */
            /*         2 == y-pixel coord in bottommost margin        */
            /* Given (ix, iy):                                        */
            /*    'startblk[ix][iy]' == starting neighbor index (sbi) */
            /*    'endblk[ix][iy]'   == ending neighbor index (ebi)   */
            /*    so that neighbors begin to be analized from index   */
            /*    'sbi' to 'ebi'.                                     */
            /* Ex. (ix, iy) = (2, 0)                                  */
            /*    ix==2 ==> x-pixel coord in rightmost margin         */
            /*    iy==0 ==> y-pixel coord in topmost margin           */
            /*    X - marks the region in the current block           */
            /*        corresponding to (ix=2, iy=0).                  */
            /*    sbi = 0 = startblk[2][0]                            */
            /*    ebi = 2 = endblk[2][0]                              */
            /*    so neighbors are analized on index range [0..2]     */
            /*                                |                       */
            /*                 nbr block 0    |  nbr block 1          */
            /*      --------------------------+------------           */
            /*           top margin      | X  |                       */
            /*      _._._._._._._._._._._._._.|                       */
            /*                           |    |                       */
            /*          current block    .r  m|  nbr block 2          */
            /*                           |i  a|                       */
            /*                           .g  g|                       */
            /*                           |h  i|                       */
            /*                           .t  n|                       */
            /*                           |    |                       */

            /* LUT for starting neighbor index given (ix, iy).        */
            static int startblk[9] = { 6, 0, 0,
                                        6,-1, 2,
                                        4, 4, 2 };
            /* LUT for ending neighbor index given (ix, iy).          */
            static int endblk[9] =   { 8, 0, 2,
                                        6,-1, 2,
                                        6, 4, 4 };

            /* Pixel coord offsets specifying the order in which neighboring */
            /* blocks are searched.  The current block is in the middle of   */
            /* 8 surrounding neighbors.  The following illustrates the order */
            /* of neighbor indices.  (Note that 9 overlaps 1.)               */
            /*                        8                                      */
            /*                      7 0 1                                    */
            /*                      6 C 2                                    */
            /*                      5 4 3                                    */
            /*                                                               */
            /*                       0  1  2  3  4  5  6  7  8                    */
            static int blkdx[9] = {  0, 1, 1, 1, 0,-1,-1,-1, 0 };  /* Delta-X     */
            static int blkdy[9] = { -1,-1, 0, 1, 1, 1, 0,-1,-1 };  /* Delta-Y     */

            /* If the margin covers more than the entire block ... */
            if (lqmParams.inv_block_margin > (lqmParams.blocksize>>1)) {
                /* Then treat this as an error. */
                throw std::domain_error("Invalid LQMParams: margin too large for blocksize (inv_block_margin > blocksize");
            }

            /* Compute the low and high pixel margin boundaries (ex. 6 pixels wide) */
            /* in the block.                                                        */
            int lo_margin = lqmParams.inv_block_margin;
            int hi_margin = lqmParams.blocksize - lqmParams.inv_block_margin - 1;

            int i = 0;
            /* Foreach minutia remaining in the list ... */
            while (i < static_cast<int>(minutiae.size())) {
                /* Assign temporary minutia pointer. */
                const Minutia& minutia = minutiae[static_cast<std::size_t>(i)];

                /* Compute block coords from minutia's pixel location. */
                int bx = minutia.x / lqmParams.blocksize;
                int by = minutia.y / lqmParams.blocksize;

                /* Compute pixel offset into the image block corresponding to the */
                /* minutia's pixel location.                                      */
                /* NOTE: The margins used here will not necessarily correspond to */
                /* the actual block boundaries used to compute the map values.    */
                /* This will be true when the image width and/or height is not an */
                /* even multiple of 'blocksize' and we are processing minutia     */
                /* located in the right-most column (or bottom-most row) of       */
                /* blocks.  I don't think this will pose a problem in practice.   */
                int px = minutia.x % lqmParams.blocksize;
                int py = minutia.y % lqmParams.blocksize;

                /* Determine if x pixel offset into the block is in the margins. */
                int ix;
                /* If x pixel offset is in left margin ... */
                if (px < lo_margin) {
                    ix = 0;
                }
                /* If x pixel offset is in right margin ... */
                else if (px > hi_margin) {
                    ix = 2;
                }
                /* Otherwise, x pixel offset is in middle of block. */
                else {
                    ix = 1;
                }

                /* Determine if y pixel offset into the block is in the margins. */
                int iy;
                /* If y pixel offset is in top margin ... */
                if (py < lo_margin) {
                    iy = 0;
                }
                /* If y pixel offset is in bottom margin ... */
                else if (py > hi_margin) {
                    iy = 2;
                }
                /* Otherwise, y pixel offset is in middle of block. */
                else {
                    iy = 1;
                }

                /* Set remove flag to false. */
                bool removed = false;

                /* If one of the minutia's pixel offsets is in a margin ... */
                if ((ix != 1) || (iy != 1)) {

                    /* Compute the starting neighbor block index for processing. */
                    int sbi = *(startblk+(iy*3)+ix);
                    /* Compute the ending neighbor block index for processing. */
                    int ebi = *(endblk+(iy*3)+ix);

                    /* Foreach neighbor in the range to be processed ... */
                    for (int ni = sbi; ni <= ebi; ++ni) {
                        /* Compute the neighbor's block coords relative to */
                        /* the block the current minutia is in.            */
                        int nbx = bx + blkdx[ni];
                        int nby = by + blkdy[ni];

                        /* If neighbor's block coords are outside of map boundaries... */
                        if ((nbx < 0) || (nbx >= mw) || (nby < 0) || (nby >= mh)) {
                            /* Then the minutia is in a margin adjacent to the edge of */
                            /* the image.                                              */
                            /* NOTE: This is true when the image width and/or height   */
                            /* is an even multiple of blocksize.  When the image is not*/
                            /* an even multiple, then some minutia may not be detected */
                            /* as being in the margin of "the image" (not the block).  */
                            /* In practice, I don't think this will impact performance.*/
                            minutiae.erase(minutiae.begin() + i);
                            /* Set remove flag to TURE. */
                            removed = true;
                            /* Break out of neighboring block loop. */
                            break;
                        }
                        /* If the neighboring block has INVALID direction ... */
                        else if (*(direction_map+(nby*mw)+nbx) == LQM_INVALID_DIR) {
                            /* Count the number of valid blocks neighboring */
                            /* the current neighbor.                        */
                            int nvalid = OpenLQM::Core::NumValid8Neigh(direction_map, nbx, nby, mw, mh);
                            /* If the number of valid neighbors is < threshold */
                            /* (ex. 7)...                                      */
                            if (nvalid < lqmParams.rm_valid_nbr_min) {
                                /* Then remove the current minutia from the list. */
                                minutiae.erase(minutiae.begin() + i);
                                /* Set remove flag to TURE. */
                                removed = true;
                                /* Break out of neighboring block loop. */
                                break;
                            }
                            /* Otherwise enough valid neighbors, so don't remove minutia */
                            /* based on this neighboring block.                          */
                        }
                        /* Otherwise neighboring block has valid direction,         */
                        /* so don't remove minutia based on this neighboring block. */
                    }

                } /* Otherwise not in margin, so skip to next minutia in list. */

                /* If current minutia not removed ... */
                if (!removed) {
                    /* Advance to the next minutia in the list. */
                    ++i;
                }
                /* Otherwise the next minutia has slid into the spot where current */
                /* minutia was removed, so don't bump minutia index.               */
            } /* End minutia loop */
        }

        void RemoveOrAdjustSideMinutiae(Minutiae& minutiae,
                        unsigned char *bdata, const int iw, const int ih,
                        int *direction_map, const int mw,
                        const LQMParams& lqmParams)
        {
            /* Allocate working memory for holding rotated y-coord of a */
            /* minutia's contour.                                       */
            int rot_y_size = (lqmParams.side_half_contour << 1) + 1;
            std::vector<int> rot_y(static_cast<std::size_t>(rot_y_size));

            /* Compute factor for converting integer directions to radians. */
            double pi_factor = LQM_PI / static_cast<double>(lqmParams.num_directions);

            int i = 0;
            /* Foreach minutia remaining in list ... */
            while (i < static_cast<int>(minutiae.size())) {
                Minutia& minutia = minutiae[static_cast<std::size_t>(i)];

                /* Extract a contour centered on the minutia point (ex. 7 pixels */
                /* in both directions).                                          */
                Contour contour;
                int ret = GetCenteredContour(contour,
                                    lqmParams.side_half_contour,
                                    minutia.x, minutia.y, minutia.ex, minutia.ey,
                                    bdata, iw, ih);

                /* If we didn't succeed in extracting a complete contour for any */
                /* other reason ...                                              */
                if ((ret == LQM_LOOP_FOUND) ||
                    (ret == LQM_IGNORE) ||
                    (ret == LQM_INCOMPLETE)) {

                    /* Remove minutia from list. */
                    minutiae.erase(minutiae.begin() + i);
                    /* No need to advance because next minutia has "slid" */
                    /* into position pointed to by 'i'.                   */
                }
                /* Otherwise, a complete contour was found and extracted ... */
                else {
                    /* Rotate contour points by negative angle of feature's direction. */
                    /* The contour of a well-formed minutia point will form a bowl     */
                    /* shape concaved in the direction of the minutia.  By rotating    */
                    /* the contour points by the negative angle of feature's direction */
                    /* the bowl will be transformed to be concaved upwards and minima  */
                    /* and maxima of the transformed y-coords can be analyzed to       */
                    /* determine if the minutia is "well-formed" or not.  If well-     */
                    /* formed then the position of the minutia point is adjusted.  If  */
                    /* not well-formed, then the minutia point is removed altogether.  */

                    /* Normal rotation of T degrees around the origin of */
                    /*      the point (x,y):                             */
                    /*         rx = x*cos(T) - y*sin(T)                  */
                    /*         ry = x*cos(T) + y*sin(T)                  */
                    /*      The rotation here is for -T degrees:         */
                    /*         rx = x*cos(-T) - y*sin(-T)                */
                    /*         ry = x*cos(-T) + y*sin(-T)                */
                    /*      which can be written:                        */
                    /*         rx = x*cos(T) + y*sin(T)                  */
                    /*         ry = x*sin(T) - y*cos(T)                  */

                    /* Convert minutia's direction to radians. */
                    double theta = static_cast<double>(minutia.direction) * pi_factor;
                    /* Compute sine and cosine values at theta for rotation. */
                    double sin_theta = std::sin(theta);
                    double cos_theta = std::cos(theta);

                    for (int j = 0; j < contour.ncontour; ++j) {
                        /* We only need to rotate the y-coord (don't worry     */
                        /* about rotating the x-coord or contour edge pixels). */
                        double drot_y = (static_cast<double>(contour.xV[static_cast<std::size_t>(j)]) * sin_theta) - (static_cast<double>(contour.yV[static_cast<std::size_t>(j)]) * cos_theta);
                        /* Need to truncate precision so that answers are consistent */
                        /* on different computer architectures when rounding doubles. */
                        drot_y = trunc_dbl_precision(drot_y);
                        rot_y[static_cast<std::size_t>(j)] = sround(drot_y);
                    }

                    /* Locate relative minima and maxima in vector of rotated */
                    /* y-coords of current minutia's contour.                 */
                    std::vector<Extremum> extrema;
                    GetExtrema(extrema, rot_y.data(), contour.ncontour);

                    /* If one and only one minima was found in rotated y-coord */
                    /* of contour ...                                          */
                    if (extrema.size() == 1 && extrema[0].type == ExtremumType::MINIMUM) {
                        /* Reset loation of minutia point to contour point at minima. */
                        minutia.x = contour.xV[static_cast<std::size_t>(extrema[0].sourceIndex)];
                        minutia.y = contour.yV[static_cast<std::size_t>(extrema[0].sourceIndex)];
                        minutia.ex = contour.exV[static_cast<std::size_t>(extrema[0].sourceIndex)];
                        minutia.ey = contour.eyV[static_cast<std::size_t>(extrema[0].sourceIndex)];

                        /* Must check if adjusted minutia is now in INVALID block ... */
                        int bx = minutia.x / lqmParams.blocksize;
                        int by = minutia.y / lqmParams.blocksize;
                        if (*(direction_map+(by*mw)+bx) == LQM_INVALID_DIR) {
                            /* Remove minutia from list. */
                            minutiae.erase(minutiae.begin() + i);
                            /* No need to advance because next minutia has "slid" */
                            /* into position pointed to by 'i'.                   */
                        }
                        else {
                            /* Advance to the next minutia in the list. */
                            ++i;
                        }

                    }
                    /* If exactly 3 min/max found and they are min-max-min ... */
                    else if (extrema.size() == 3 && extrema[0].type == ExtremumType::MINIMUM) {
                        /* Choose minima location with smallest rotated y-coord. */
                        int minloc;
                        if (extrema[0].val < extrema[2].val) {
                            minloc = extrema[0].sourceIndex;
                        }
                        else {
                            minloc = extrema[2].sourceIndex;
                        }

                        /* Reset loation of minutia point to contour point at minima. */
                        minutia.x = contour.xV[static_cast<std::size_t>(minloc)];
                        minutia.y = contour.yV[static_cast<std::size_t>(minloc)];
                        minutia.ex = contour.exV[static_cast<std::size_t>(minloc)];
                        minutia.ey = contour.eyV[static_cast<std::size_t>(minloc)];

                        /* Must check if adjusted minutia is now in INVALID block ... */
                        int bx = minutia.x / lqmParams.blocksize;
                        int by = minutia.y / lqmParams.blocksize;
                        if (*(direction_map+(by*mw)+bx) == LQM_INVALID_DIR) {
                            /* Remove minutia from list. */
                            minutiae.erase(minutiae.begin() + i);

                            /* No need to advance because next minutia has "slid" */
                            /* into position pointed to by 'i'.                   */
                        }
                        else {
                            /* Advance to the next minutia in the list. */
                            i++;
                        }
                    }
                    /* Otherwise, ... */
                    else {
                        /* Remove minutia from list. */
                        minutiae.erase(minutiae.begin() + i);

                        /* No need to advance because next minutia has "slid" */
                        /* into position pointed to by 'i'.                   */
                    }
                } /* End else contour extracted. */
            } /* End while not end of minutiae list. */
        }

        int GetCenteredContour(Contour& ocontour,
                        const int half_contour,
                        const int x_loc, const int y_loc,
                        const int x_edge, const int y_edge,
                        unsigned char *bdata, const int iw, const int ih)
        {
            /* Compute maximum length of complete contour */
            /* (2 half contours + feature point).         */
            int max_contour = (half_contour << 1) + 1;

            /* Get 1st half contour with clockwise neighbor trace. */
            Contour half1;
            int ret = TraceContour(half1,
                            half_contour, x_loc, y_loc, x_loc, y_loc, x_edge, y_edge,
                            LQM_SCAN_CLOCKWISE, bdata, iw, ih);

            /* If trace was not possible ... */
            if (ret == LQM_IGNORE) {
                /* Return LQM_IGNORE */
                return LQM_IGNORE;
            }

            /* If 1st half contour forms a loop ... */
            if (ret == LQM_LOOP_FOUND) {
                /* Return LQM_LOOP_FOUND */
                return LQM_LOOP_FOUND;
            }

            /* If 1st half contour not complete ... */
            if (half1.ncontour < half_contour) {
                /* Return with contour length equal to 0. */
                return LQM_INCOMPLETE;
            }

            /* Otherwise, we have a complete 1st half contour...           */
            /* Get 2nd half contour with counter-clockwise neighbor trace. */
            /* Use the last point from the first contour trace as the      */
            /* point to test for a loop when tracing the second contour.   */
            Contour half2;
            ret = TraceContour(half2,
                            half_contour, half1.xV[static_cast<std::size_t>(half1.ncontour - 1)], half1.yV[static_cast<std::size_t>(half1.ncontour - 1)],
                            x_loc, y_loc, x_edge, y_edge,
                            LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih);

            /* If 2nd trace was not possible ... */
            if (ret == LQM_IGNORE) {
                /* Return with contour length equal to 0. */
                return LQM_IGNORE;
            }

            /* If 2nd trace forms a loop ... */
            if (ret == LQM_LOOP_FOUND) {
                /* Return LQM_LOOP_FOUND */
                return LQM_LOOP_FOUND;
            }

            /* If 2nd half contour not complete ... */
            if (half2.ncontour < half_contour) {
                /* Return with contour length equal to 0. */
                return LQM_INCOMPLETE;
            }

            /* Otherwise we have a full 1st half contour and a 2nd half contour */
            /* that do not form a loop and are complete.  We now need to        */
            /* concatenate the two half contours into one longer contour.       */

            /* Allocate output contour list. */
            ocontour.Reserve(max_contour);

            /* Copy 1st half contour into output contour buffers.      */
            /* This contour was collected clockwise, so it's points    */
            /* are entered in reverse order of the trace.  The result  */
            /* is the first point in the output contour if farthest    */
            /* from the starting feature point.                        */
            for (int j = half1.ncontour - 1; j >= 0; --j) {
                ocontour.AddPoint(half1.xV[static_cast<std::size_t>(j)], half1.yV[static_cast<std::size_t>(j)], half1.exV[static_cast<std::size_t>(j)], half1.eyV[static_cast<std::size_t>(j)]);
            }

            /* Next, store starting feature point into output contour buffers. */
            ocontour.AddPoint(x_loc, y_loc, x_edge, y_edge);

            /* Now, append 2nd half contour to permanent contour buffers.  */
            for (int i = 0; i < half2.ncontour; ++i) {
                ocontour.AddPoint(half2.xV[static_cast<std::size_t>(i)], half2.yV[static_cast<std::size_t>(i)], half2.exV[static_cast<std::size_t>(i)], half2.eyV[static_cast<std::size_t>(i)]);
            }

            /* Return normally. */
            return 0;
        }

        void GetExtrema(std::vector<Extremum>& extrema, const int *items, const int num)
        {
            /* Determine maximum length for allocation of buffers. */

            /* If there are fewer than 3 items ...                */
            if (num < 3) {
                /* Then no min/max is possible, so return */
                return;
            }
            /* Otherwise, set allocation length to number of items - 2    */
            /* (one for the first item in the list, and on for the last). */
            /* Every other intermediate point can potentially represent a */
            /* min or max.                                                */
            int minmax_alloc = num - 2;
            /* Allocate the buffers. */
            extrema.reserve(static_cast<std::size_t>(minmax_alloc));

            /* Start witht the first item in the list. */
            int i = 0;

            /* Get starting state between first pair of items. */
            int diff = items[1] - items[0];
            ExtremumType state;
            if (diff > 0) {
                state = ExtremumType::MAXIMUM;
            }
            else if (diff < 0) {
                state = ExtremumType::MINIMUM;
            }
            else {
                state = ExtremumType::INDETERMINATE;
            }

            /* Set start location to first item in list. */
            int start = 0;

            /* Bump to next item in list. */
            ++i;

            int loc;
            /* While not at the last item in list. */
            while (i < num - 1) {
                /* Compute difference between next pair of items. */
                diff = items[i+1] - items[i];
                /* If items are increasing ... */
                if (diff > 0) {
                    /* If previously increasing ... */
                    if (state == ExtremumType::MAXIMUM) {
                        /* Reset start to current location. */
                        start = i;
                    }
                    /* If previously decreasing ... */
                    else if (state == ExtremumType::MINIMUM) {
                        /* Then we have incurred a minima ... */
                        /* Compute midpoint of minima. */
                        loc = (start + i)/2;
                        /* Store minima data */
                        extrema.emplace_back(items[loc], ExtremumType::MINIMUM, loc);
                        /* Change state to increasing. */
                        state = ExtremumType::MAXIMUM;
                        /* Reset start location. */
                        start = i;
                    }
                    /* If previously level (this state only can occur at the */
                    /* beginning of the list of items) ...                   */
                    else {
                        /* If more than one level state in a row ... */
                        if (i - start > 1) {
                            /* Then consider a minima ... */
                            /* Compute midpoint of minima. */
                            loc = (start + i)/2;
                            /* Store minima data. */
                            extrema.emplace_back(items[loc], ExtremumType::MINIMUM, loc);
                            /* Change state to increasing. */
                            state = ExtremumType::MAXIMUM;
                            /* Reset start location. */
                            start = i;
                        }
                        /* Otherwise, ignore single level state. */
                        else {
                            /* Change state to increasing. */
                            state = ExtremumType::MAXIMUM;
                            /* Reset start location. */
                            start = i;
                        }
                    }
                }
                /* If items are decreasing ... */
                else if (diff < 0) {
                    /* If previously decreasing ... */
                    if (state == ExtremumType::MINIMUM) {
                        /* Reset start to current location. */
                        start = i;
                    }
                    /* If previously increasing ... */
                    else if (state == ExtremumType::MAXIMUM) {
                        /* Then we have incurred a maxima ... */
                        /* Compute midpoint of maxima. */
                        loc = (start + i)/2;
                        /* Store maxima data. */
                        extrema.emplace_back(items[loc], ExtremumType::MAXIMUM, loc);
                        /* Change state to decreasing. */
                        state = ExtremumType::MINIMUM;
                        /* Reset start location. */
                        start = i;
                    }
                    /* If previously level (this state only can occur at the */
                    /* beginning of the list of items) ...                   */
                    else {
                        /* If more than one level state in a row ... */
                        if (i - start > 1) {
                            /* Then consider a maxima ... */
                            /* Compute midpoint of maxima. */
                            loc = (start + i)/2;
                            /* Store maxima data. */
                            extrema.emplace_back(items[loc], ExtremumType::MAXIMUM, loc);
                            /* Change state to decreasing. */
                            state = ExtremumType::MINIMUM;
                            /* Reset start location. */
                            start = i;
                        }
                        /* Otherwise, ignore single level state. */
                        else {
                            /* Change state to decreasing. */
                            state = ExtremumType::MINIMUM;
                            /* Reset start location. */
                            start = i;
                        }
                    }
                }
                /* Otherwise, items are level, so continue to next item pair. */

                /* Advance to next item pair in list. */
                i++;
            }
        }

        void RemoveHooks(Minutiae& minutiae,
                        unsigned char *bdata, const int iw, const int ih,
                        const LQMParams& lqmParams)
        {
            /* Allocate list of minutia indices that upon completion of testing */
            /* should be removed from the minutiae lists. Initialized to false  */
            std::vector<bool> to_remove(minutiae.size(), false);

            /* Compute number directions in full circle. */
            int full_ndirs = lqmParams.num_directions << 1;
            /* Compute number of directions in 45=(180/4) degrees. */
            int qtr_ndirs = lqmParams.num_directions >> 2;

            /* Minimum allowable deltadir to consider joining minutia.               */
            /* (The closer the deltadir is to 180 degrees, the more likely the join. */
            /* When ndirs==16, then this value is 11=(3*4)-1 == 123.75 degrees.      */
            /* I chose to parameterize this threshold based on a fixed fraction of   */
            /* 'ndirs' rather than on passing in a parameter in degrees and doing    */
            /* the conversion.  I doubt the difference matters.                      */
            int min_deltadir = (3 * qtr_ndirs) - 1;

            int f = 0;
            /* Foreach primary (first) minutia (except for last one in list) ... */
            while (f < (static_cast<int>(minutiae.size())) - 1) {

                /* If current first minutia not previously set to be removed. */
                if (!to_remove[static_cast<std::size_t>(f)]) {
                    /* Set first minutia to temporary pointer. */
                    Minutia& minutia1 = minutiae[static_cast<std::size_t>(f)];
                    /* Foreach secondary (second) minutia to right of first minutia ... */
                    int s = f + 1;
                    while (s < static_cast<int>(minutiae.size())) {
                        /* Set second minutia to temporary pointer. */
                        Minutia& minutia2 = minutiae[static_cast<std::size_t>(s)];

                        /* The binary image is potentially being edited during each */
                        /* iteration of the secondary minutia loop, therefore       */
                        /* minutia pixel values may be changed.  We need to catch   */
                        /* these events by using the next 2 tests.                  */

                        /* If the first minutia's pixel has been previously changed... */
                        if (*(bdata+(minutia1.y*iw)+minutia1.x) != minutia1.type) {
                            /* Then break out of secondary loop and skip to next first. */
                            break;
                        }

                        /* If the second minutia's pixel has been previously changed... */
                        if (*(bdata+(minutia2.y*iw)+minutia2.x) != minutia2.type) {
                            /* Set to remove second minutia. */
                            to_remove[static_cast<std::size_t>(s)] = true;
                        }

                        /* If the second minutia not previously set to be removed. */
                        if (!to_remove[static_cast<std::size_t>(s)]) {
                            /* Compute delta y between 1st & 2nd minutiae and test. */
                            int delta_y = minutia2.y - minutia1.y;
                            /* If delta y small enough (ex. < 8 pixels) ... */
                            if (delta_y <= lqmParams.max_rmtest_dist) {
                                /* Compute Euclidean distance between 1st & 2nd mintuae. */
                                double dist = std::hypot(minutia1.x - minutia2.x, minutia1.y - minutia2.y);
                                /* If distance is NOT too large (ex. < 8 pixels) ... */
                                if (dist <= lqmParams.max_rmtest_dist) {
                                    /* Compute "inner" difference between directions on */
                                    /* a full circle and test.                          */
                                    int deltadir = OpenLQM::Core::ClosestDirDist(minutia1.direction, minutia2.direction, full_ndirs);
                                    if (deltadir == LQM_INVALID_DIR) {
                                        throw std::invalid_argument("RemoveHooks was passed a minutia with invalid direction");
                                    }
                                    /* If the difference between dirs is large enough ...  */
                                    /* (the more 1st & 2nd point away from each other the  */
                                    /* more likely they should be joined)                  */
                                    if (deltadir > min_deltadir) {
                                        /* If 1st & 2nd minutiae are NOT same type ... */
                                        if (minutia1.type != minutia2.type) {
                                            /* Check to see if pair on a hook with contour */
                                            /* of specified length (ex. 15 pixels) ...     */

                                            int ret = OnHook(minutia1, minutia2,
                                                            lqmParams.max_hook_len,
                                                            bdata, iw, ih);

                                            /* If hook detected between pair ... */
                                            if (ret == LQM_HOOK_FOUND) {
                                                /* Set to remove first minutia. */
                                                to_remove[static_cast<std::size_t>(f)] = true;
                                                /* Set to remove second minutia. */
                                                to_remove[static_cast<std::size_t>(s)] = true;
                                            }
                                            /* If hook test IGNORED ... */
                                            else if (ret == LQM_IGNORE) {
                                                /* Set to remove first minutia. */
                                                to_remove[static_cast<std::size_t>(f)] = true;
                                                /* Skip to next 1st minutia by breaking out of */
                                                /* inner secondary loop.                       */
                                                break;
                                            }
                                            /* Otherwise, no hook found, so skip to next */
                                            /* second minutia.                           */
                                        }
                                        /* End different type test. */
                                    }/* End deltadir test. */
                                }/* End distance test. */
                            }
                            /* Otherwise, current 2nd too far below 1st, so skip to next */
                            /* 1st minutia.                                              */
                            else {
                                /* Break out of inner secondary loop. */
                                break;
                            }/* End delta-y test. */
                        }/* End if !to_remove[s] */

                        /* Bump to next second minutia in minutiae list. */
                        ++s;
                    }/* End secondary minutiae loop. */

                }/* Otherwise, first minutia already flagged to be removed. */

                /* Bump to next first minutia in minutiae list. */
                ++f;
            }/* End primary minutiae loop. */

            /* Now remove all minutiae in list that have been flagged for removal. */
            /* NOTE: Need to remove the minutia from their lists in reverse       */
            /*       order, otherwise, indices will be off.                       */
            for (int i = (static_cast<int>(minutiae.size())) - 1; i >= 0; --i) {
                /* If the current minutia index is flagged for removal ... */
                if (to_remove[static_cast<std::size_t>(i)]) {
                    /* Remove the minutia from the minutiae list. */
                    minutiae.erase(minutiae.begin() + i);
                }
            }
        }

        int OnHook(const Minutia& minutia1, const Minutia& minutia2,
                    const int max_hook_len,
                    unsigned char *bdata, const int iw, const int ih)
        {
            /* NOTE: This routine should only be called when the 2 minutia points */
            /*       are of "opposite" type.                                      */

            /* Trace the contour of the feature starting at the 1st minutia's     */
            /* "edge" point and stepping along up to the specified maximum number */
            /* of steps or until the 2nd minutia point is encountered.            */
            /* First search for edge neighbors clockwise.                         */
            int ret;
            {
                Contour contour;
                ret = TraceContour(contour, max_hook_len,
                                    minutia2.x, minutia2.y, minutia1.ex, minutia1.ey,
                                    minutia1.x, minutia1.y,
                                    LQM_SCAN_CLOCKWISE, bdata, iw, ih);
            }

            /* If trace was not possible, return LQM_IGNORE. */
            if (ret == LQM_IGNORE) {
                return ret;
            }

            /* If the trace encountered the second minutia point ... */
            if (ret == LQM_LOOP_FOUND) {
                return LQM_HOOK_FOUND;
            }

            /* Otherwise, the trace successfully followed the contour, but did */
            /* not encounter the 2nd minutia point within the specified number */
            /* of steps.                                                       */

            /* Try searching contour from 1st minutia "edge" searching for */
            /* edge neighbors counter-clockwise.                           */
            {
                Contour contour;
                ret = TraceContour(contour, max_hook_len,
                                    minutia2.x, minutia2.y, minutia1.ex, minutia1.ey,
                                    minutia1.x, minutia1.y,
                                    LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih);
            }
            /* If trace was not possible, return LQM_IGNORE. */
            if (ret == LQM_IGNORE) {
                return ret;
            }

            /* If the trace encountered the second minutia point ... */
            if (ret == LQM_LOOP_FOUND) {
                return LQM_HOOK_FOUND;
            }

            /* The trace successfully followed the 1st minutia's contour, but   */
            /* did not encounter the 2nd minutia point within the specified number */
            /* of steps ...      */
            /* Return hook NOT found (LQM_FALSE). */
            return LQM_FALSE;
        }

        void RemoveOverlaps(Minutiae& minutiae,
                            unsigned char *bdata, const int iw,
                            const LQMParams& lqmParams)
        {
            /* Allocate list of minutia indices that upon completion of testing */
            /* should be removed from the minutiae lists. Initialized to false. */
            std::vector<bool> to_remove(minutiae.size(), false);

            /* Compute number directions in full circle. */
            int full_ndirs = lqmParams.num_directions<<1;
            /* Compute number of directions in 45=(180/4) degrees. */
            int qtr_ndirs = lqmParams.num_directions>>2;
            /* Compute number of directions in 90=(180/2) degrees. */
            int half_ndirs = lqmParams.num_directions>>1;

            /* Minimum allowable deltadir to consider joining minutia.               */
            /* (The closer the deltadir is to 180 degrees, the more likely the join. */
            /* When ndirs==16, then this value is 11=(3*4)-1 == 123.75 degrees.      */
            /* I chose to parameterize this threshold based on a fixed fraction of   */
            /* 'ndirs' rather than on passing in a parameter in degrees and doing    */
            /* the conversion.  I doubt the difference matters.                      */
            int min_deltadir = (3 * qtr_ndirs) - 1;

            int f = 0;
            /* Foreach primary (first) minutia (except for last one in list) ... */
            while (f < (static_cast<int>(minutiae.size())) - 1) {
                /* If current first minutia not previously set to be removed. */
                if (!to_remove[static_cast<std::size_t>(f)]) {
                    /* Set first minutia to temporary pointer. */
                    Minutia& minutia1 = minutiae[static_cast<std::size_t>(f)];
                    /* Foreach secondary (second) minutia to right of first minutia ... */
                    int s = f+1;
                    while (s < static_cast<int>(minutiae.size())) {
                        /* Set second minutia to temporary pointer. */
                        Minutia& minutia2 = minutiae[static_cast<std::size_t>(s)];

                        /* The binary image is potentially being edited during each */
                        /* iteration of the secondary minutia loop, therefore       */
                        /* minutia pixel values may be changed.  We need to catch   */
                        /* these events by using the next 2 tests.                  */

                        /* If the first minutia's pixel has been previously changed... */
                        if (*(bdata+(minutia1.y*iw)+minutia1.x) != minutia1.type) {
                            /* Then break out of secondary loop and skip to next first. */
                            break;
                        }

                        /* If the second minutia's pixel has been previously changed... */
                        if (*(bdata+(minutia2.y*iw)+minutia2.x) != minutia2.type) {
                            /* Set to remove second minutia. */
                            to_remove[static_cast<std::size_t>(s)] = true;
                        }

                        /* If the second minutia not previously set to be removed. */
                        if (!to_remove[static_cast<std::size_t>(s)]) {
                            /* Compute delta y between 1st & 2nd minutiae and test. */
                            int delta_y = minutia2.y - minutia1.y;
                            /* If delta y small enough (ex. < 8 pixels) ... */
                            if (delta_y <= lqmParams.max_overlap_dist) {
                                /* Compute Euclidean distance between 1st & 2nd mintuae. */
                                double dist = std::hypot(minutia1.x - minutia2.x, minutia1.y - minutia2.y);

                                /* If distance is NOT too large (ex. < 8 pixels) ... */
                                if (dist <= lqmParams.max_overlap_dist) {
                                    /* Compute "inner" difference between directions on */
                                    /* a full circle and test.                          */
                                    int deltadir = OpenLQM::Core::ClosestDirDist(minutia1.direction, minutia2.direction, full_ndirs);
                                    if (deltadir == LQM_INVALID_DIR) {
                                        throw std::invalid_argument("RemoveOverlaps was passed a minutia with invalid direction");
                                    }
                                    /* If the difference between dirs is large enough ...  */
                                    /* (the more 1st & 2nd point away from each other the  */
                                    /* more likely they should be joined)                  */
                                    if (deltadir > min_deltadir) {
                                        /* If 1st & 2nd minutiae are same type ... */
                                        if (minutia1.type == minutia2.type) {
                                            /* Test to see if both are on opposite sides */
                                            /* of an overlap.                            */

                                            /* Compute direction of "joining" vector.      */
                                            /* First, compute direction of line from first */
                                            /* to second minutia points.                   */
                                            int joindir = LineToDirection(minutia1.x, minutia1.y,
                                                                        minutia2.x, minutia2.y,
                                                                        lqmParams.num_directions);

                                            /* Comptue opposite direction of first minutia. */
                                            int opp1dir = (minutia1.direction + lqmParams.num_directions) % full_ndirs;
                                            /* Take "inner" distance on full circle between */
                                            /* the first minutia's opposite direction and   */
                                            /* the joining direction.                       */
                                            joindir = std::abs(opp1dir - joindir);
                                            joindir = std::min(joindir, full_ndirs - joindir);

                                            /* If the joining angle is <= 90 degrees OR   */
                                            /*    the 2 points are sufficiently close AND */
                                            /*    a free path exists between pair ...     */
                                            if (
                                                joindir <= half_ndirs
                                                || (dist <= lqmParams.max_overlap_join_dist
                                                && OpenLQM::Core::FreePath(
                                                    minutia1.x, minutia1.y, minutia2.x, minutia2.y,
                                                    bdata, iw, lqmParams
                                                ))
                                            ) {
                                                /* Then assume overlap, so ...             */
                                                /* Set to remove first minutia. */
                                                to_remove[static_cast<std::size_t>(f)] = true;
                                                /* Set to remove second minutia. */
                                                to_remove[static_cast<std::size_t>(s)] = true;
                                            }
                                            /* Otherwise, pair not on an overlap, so skip */
                                            /* to next second minutia.                    */
                                        }
                                        /* End same type test. */
                                    }/* End deltadir test. */
                                }/* End distance test. */
                            }
                            /* Otherwise, current 2nd too far below 1st, so skip to next */
                            /* 1st minutia.                                              */
                            else {
                                /* Break out of inner secondary loop. */
                                break;
                            }/* End delta-y test. */
                        }/* End if !to_remove[s] */

                        /* Bump to next second minutia in minutiae list. */
                        s++;
                    }/* End secondary minutiae loop. */

                }

                /* Bump to next first minutia in minutiae list. */
                f++;
            }/* End primary minutiae loop. */

            /* Now remove all minutiae in list that have been flagged for removal. */
            /* NOTE: Need to remove the minutia from their lists in reverse       */
            /*       order, otherwise, indices will be off.                       */
            for (int i = (static_cast<int>(minutiae.size())) - 1; i >= 0; --i) {
                /* If the current minutia index is flagged for removal ... */
                if (to_remove[static_cast<std::size_t>(i)]) {
                    /* Remove the minutia from the minutiae list. */
                    minutiae.erase(minutiae.begin() + i);
                }
            }
        }

        void RemoveMalformations(Minutiae& minutiae,
                                unsigned char *bdata, const int iw, const int ih,
                                int *low_flow_map, const int mw,
                                const LQMParams& lqmParams)
        {
            for (int i = (static_cast<int>(minutiae.size())) - 1; i >= 0; --i) {
                int ret;
                const Minutia& minutia = minutiae[static_cast<std::size_t>(i)];
                int ax1, ay1, bx1, by1;
                {
                    Contour contour;
                    ret = TraceContour(contour,
                                        lqmParams.malformation_steps_2,
                                        minutia.x, minutia.y,
                                        minutia.x, minutia.y, minutia.ex, minutia.ey,
                                        LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih);

                    /* If trace was not possible OR loop found OR */
                    /* contour is incomplete ...                  */
                    if (ret == LQM_IGNORE || ret == LQM_LOOP_FOUND || contour.ncontour < lqmParams.malformation_steps_2) {
                        /* Then remove the minutia. */
                        minutiae.erase(minutiae.begin() + i);
                        continue;
                    }
                    /* Otherwise, traced contour is complete. */

                    /* Store 'A1' contour point. */
                    ax1 = contour.xV[static_cast<std::size_t>(lqmParams.malformation_steps_1 - 1)];
                    ay1 = contour.yV[static_cast<std::size_t>(lqmParams.malformation_steps_1 - 1)];

                    /* Store 'B1' contour point. */
                    bx1 = contour.xV[static_cast<std::size_t>(lqmParams.malformation_steps_2 - 1)];
                    by1 = contour.yV[static_cast<std::size_t>(lqmParams.malformation_steps_2 - 1)];
                }

                Contour contour;
                ret = TraceContour(contour,
                                lqmParams.malformation_steps_2,
                                minutia.x, minutia.y,
                                minutia.x, minutia.y, minutia.ex, minutia.ey,
                                LQM_SCAN_CLOCKWISE, bdata, iw, ih);

                /* If trace was not possible OR loop found OR */
                /* contour is incomplete ...                  */
                if (ret == LQM_IGNORE || ret == LQM_LOOP_FOUND || contour.ncontour < lqmParams.malformation_steps_2) {
                    /* Then remove the minutia. */
                    minutiae.erase(minutiae.begin() + i);
                    continue;
                }
                /* Otherwise, traced contour is complete. */

                /* Store 'A2' contour point. */
                int ax2 = contour.xV[static_cast<std::size_t>(lqmParams.malformation_steps_1 - 1)];
                int ay2 = contour.yV[static_cast<std::size_t>(lqmParams.malformation_steps_1 - 1)];

                /* Store 'B2' contour point. */
                int bx2 = contour.xV[static_cast<std::size_t>(lqmParams.malformation_steps_2 - 1)];
                int by2 = contour.yV[static_cast<std::size_t>(lqmParams.malformation_steps_2 - 1)];

                /* Compute distances along A & B paths. */
                double a_dist = std::hypot(ax1 - ax2, ay1 - ay2);
                double b_dist = std::hypot(bx1 - bx2, by1 - by2);

                /* Compute block coords from minutia's pixel location. */
                int blk_x = minutia.x / lqmParams.blocksize;
                int blk_y = minutia.y / lqmParams.blocksize;

                /* Check to see if distances are not zero. */
                if (a_dist == 0.0 || b_dist == 0.0) {
                    /* Remove the malformation minutia. */
                    minutiae.erase(minutiae.begin() + i);
                    continue;
                }

                /* Determine if minutia is in LOW RIDGE FLOW block. */
                int fmapval = *(low_flow_map+(blk_y*mw)+blk_x);
                if (fmapval) {
                    /* If in LOW RIDGE LFOW, conduct a cursory distance test. */
                    /* Need to test this out!                                 */
                    if (b_dist > lqmParams.max_malformation_dist) {
                        /* Remove the malformation minutia. */
                        minutiae.erase(minutiae.begin() + i);
                        continue;
                    }
                }

                /* Compute points on line between the points A & B. */
                std::vector<std::pair<int, int>> points;
                OpenLQM::Core::LinePoints(points, bx1, by1, bx2, by2);

                /* Foreach remaining point along line segment ... */
                for (int j = 0; j < static_cast<int>(points.size()); ++j) {
                    /* If B path contains pixel opposite minutia type ... */
                    if (*(bdata + (points[static_cast<std::size_t>(j)].second * iw) + points[static_cast<std::size_t>(j)].first) != minutia.type) {
                        /* Compute ratio of A & B path lengths. */
                        double ratio = b_dist / a_dist;
                        /* Need to truncate precision so that answers are  */
                        /* consistent on different computer architectures. */
                        ratio = trunc_dbl_precision(ratio);
                        /* If the B path is sufficiently longer than A path ... */
                        if (ratio > lqmParams.min_malformation_ratio) {
                            /* Remove the malformation minutia. */
                            /* Then remove the minutia. */
                            minutiae.erase(minutiae.begin() + i);
                            break;
                        }
                    }
                }
            }
        }

        void RemovePores(Minutiae& minutiae,
                            unsigned char *bdata, const int iw, const int ih,
                            int *direction_map, int *low_flow_map,
                            int *high_curve_map, const int mw,
                            const LQMParams& lqmParams)
        {
            /*      This routine attempts to locate the following points on all */
            /*      minutia within the feature list.                            */
            /*      1. Compute R 3 pixels opposite the feature direction from   */
            /*         feature point F.                                         */
            /*      2. Find white pixel transitions P & Q within 12 steps from  */
            /*         from R perpendicular to the feature's direction.         */
            /*      3. Find points B & D by walking white edge from P.          */
            /*      4. Find points A & C by walking white edge from Q.          */
            /*      5. Measure squared distances between A-B and C-D.           */
            /*      6. Compute ratio of squared distances and compare against   */
            /*         threshold (2.25).  If A-B sufficiently larger than C-D,  */
            /*         then assume NOT pore, otherwise flag the feature point F.*/
            /*      If along the way, finding any of these points fails, then   */
            /*      assume the feature is a pore and flag it.                   */
            /*                                                                  */
            /*                      A                                           */
            /*                _____._                                           */
            /*                       ----___     Q      C                       */
            /*                ------____    ---_.________.___                   */
            /*                          ---_                                    */
            /*                (valley)    F.\   .R  (ridge)                     */
            /*                          ____/                                   */
            /*                ______----    ___-.--------.---                   */
            /*                       ____---     P      D                       */
            /*                -----.-                                           */
            /*                      B                                           */
            /*                                                                  */
            /*              AB^2/CD^2 <= 2.25  then flag feature                */
            /*                                                                  */

            /* Factor for converting integer directions into radians. */
            double pi_factor = LQM_PI/static_cast<double>(lqmParams.num_directions);

            /* Initialize to the beginning of the minutia list. */
            int i = 0;
            /* Foreach minutia remaining in the list ... */
            while (i < static_cast<int>(minutiae.size())) {
                /* Set temporary minutia pointer. */
                const Minutia& minutia = minutiae[static_cast<std::size_t>(i)];

                /* Compute block coords from minutia point. */
                int blk_x = minutia.x / lqmParams.blocksize;
                int blk_y = minutia.y / lqmParams.blocksize;

                /* If minutia in LOW RIDGE FLOW or HIGH CURVATURE block */
                /* with a valid direction ...                           */
                if ((*(low_flow_map+(blk_y*mw)+blk_x) ||
                    *(high_curve_map+(blk_y*mw)+blk_x)) &&
                    (*(direction_map+(blk_y*mw)+blk_x) >= 0)
                ) {
                    /* Compute radian angle from minutia direction. */
                    double theta = static_cast<double>(minutia.direction) * pi_factor;
                    /* Compute sine and cosine factors of this angle. */
                    double sin_theta = std::sin(theta);
                    double cos_theta = std::cos(theta);
                    /* Translate the minutia point (ex. 3 pixels) in opposite */
                    /* direction minutia is pointing.  Call this point 'R'.   */
                    double drx = static_cast<double>(minutia.x) -
                                (sin_theta * static_cast<double>(lqmParams.pores_trans_r));
                    double dry = static_cast<double>(minutia.y) +
                                (cos_theta * static_cast<double>(lqmParams.pores_trans_r));
                    /* Need to truncate precision so that answers are consistent */
                    /* on different computer architectures when rounding doubles. */
                    drx = trunc_dbl_precision(drx);
                    dry = trunc_dbl_precision(dry);
                    int rx = sround(drx);
                    int ry = sround(dry);

                    /* If 'R' is opposite color from minutia type ... */
                    if (*(bdata+(ry*iw)+rx) != minutia.type) {
                        /* Search a specified number of steps (ex. 12) from 'R' in a */
                        /* perpendicular direction from the minutia direction until  */
                        /* the first white pixel is found.  If a white pixel is      */
                        /* found within the specified number of steps, then call     */
                        /* this point 'P' (storing the point's edge pixel as well).  */
                        FeaturePoint p;
                        if (!SearchInDirection(p,
                                            minutia.type,
                                            rx, ry, -cos_theta, -sin_theta,
                                            lqmParams.pores_perp_steps,
                                            bdata, iw, ih))
                        {
                            /* P not found. Remove the minutia. */
                            minutiae.erase(minutiae.begin() + i);
                            continue;
                        }

                        int ret;

                        /* Trace contour from P's edge pixel in counter-clockwise  */
                        /* scan and step along specified number of steps (ex. 10). */
                        /* Store last point in contour as point 'B'.               */
                        int bx, by;
                        {
                            Contour contour;
                            ret = TraceContour(contour,
                                                lqmParams.pores_steps_fwd,
                                                p.x, p.y, p.x, p.y, p.ex, p.ey,
                                                LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih);

                            /* If trace was not possible OR loop found OR */
                            /* contour is incomplete ...                  */
                            if (ret == LQM_IGNORE || ret == LQM_LOOP_FOUND || contour.ncontour < lqmParams.pores_steps_fwd) {
                                /* Then remove the minutia. */
                                minutiae.erase(minutiae.begin() + i);
                                continue;
                            }
                            /* Otherwise, traced contour is complete. */

                            bx = contour.xV[static_cast<std::size_t>(contour.ncontour - 1)];
                            by = contour.yV[static_cast<std::size_t>(contour.ncontour - 1)];
                        }

                        /* Trace contour from P's edge pixel in clockwise scan */
                        /* and step along specified number of steps (ex. 8).   */
                        /* Store last point in contour as point 'D'.           */
                        int dx, dy;
                        {
                            Contour contour;
                            ret = TraceContour(contour,
                                                lqmParams.pores_steps_bwd,
                                                p.x, p.y, p.x, p.y, p.ex, p.ey,
                                                LQM_SCAN_CLOCKWISE, bdata, iw, ih);

                            /* If trace was not possible OR loop found OR */
                            /* contour is incomplete ...                  */
                            if (ret == LQM_IGNORE || ret == LQM_LOOP_FOUND || contour.ncontour < lqmParams.pores_steps_bwd) {
                                /* Then remove the minutia. */
                                minutiae.erase(minutiae.begin() + i);
                                continue;
                            }
                            /* Otherwise, traced contour is complete. */

                            dx = contour.xV[static_cast<std::size_t>(contour.ncontour - 1)];
                            dy = contour.yV[static_cast<std::size_t>(contour.ncontour - 1)];
                        }

                        /* Search a specified number of steps (ex. 12) from */
                        /* 'R' in opposite direction of that used to find   */
                        /* 'P' until the first white pixel is found.  If a  */
                        /* white pixel is found within the specified number */
                        /* of steps, then call this point 'Q' (storing the  */
                        /* point's edge pixel as well).                     */
                        FeaturePoint q;
                        if (!SearchInDirection(q,
                                                minutia.type,
                                                rx, ry, cos_theta, sin_theta,
                                                lqmParams.pores_perp_steps,
                                                bdata, iw, ih))
                        {
                            /* Q not found. Remove the minutia. */
                            minutiae.erase(minutiae.begin() + i);
                            continue;
                        }

                        /* Trace contour from Q's edge pixel in clockwise */
                        /* scan and step along specified number of steps  */
                        /* (ex. 10).                                      */
                        /* Store last point in contour as point 'A'.      */
                        int ax, ay;
                        {
                            Contour contour;
                            ret = TraceContour(contour,
                                        lqmParams.pores_steps_fwd,
                                        q.x, q.y, q.x, q.y, q.ex, q.ey,
                                        LQM_SCAN_CLOCKWISE, bdata, iw, ih);

                            /* If trace was not possible OR loop found OR */
                            /* contour is incomplete ...                  */
                            if (ret == LQM_IGNORE || ret == LQM_LOOP_FOUND || contour.ncontour < lqmParams.pores_steps_fwd) {
                                /* Then remove the minutia. */
                                minutiae.erase(minutiae.begin() + i);
                                continue;
                            }
                            /* Otherwise, traced contour is complete. */

                            ax = contour.xV[static_cast<std::size_t>(contour.ncontour - 1)];
                            ay = contour.yV[static_cast<std::size_t>(contour.ncontour - 1)];
                        }

                        /* Trace contour from Q's edge pixel in    */
                        /* counter-clockwise scan and step along a */
                        /* specified number of steps (ex. 8).      */
                        /* Store last point in contour as 'C'.     */
                        int cx, cy;
                        {
                            Contour contour;
                            ret = TraceContour(contour,
                                        lqmParams.pores_steps_bwd,
                                        q.x, q.y, q.x, q.y, q.ex, q.ey,
                                        LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih);

                            /* If trace was not possible OR loop found OR */
                            /* contour is incomplete ...                  */
                            if (ret == LQM_IGNORE || ret == LQM_LOOP_FOUND || contour.ncontour < lqmParams.pores_steps_bwd) {
                                /* Then remove the minutia. */
                                minutiae.erase(minutiae.begin() + i);
                                continue;
                            }
                            /* Otherwise, traced contour is complete. */

                            cx = contour.xV[static_cast<std::size_t>(contour.ncontour - 1)];
                            cy = contour.yV[static_cast<std::size_t>(contour.ncontour - 1)];
                        }

                        /* Compute squared distance between points */
                        /* 'A' and 'B'.                            */
                        double ab2 = SquaredDistance(ax, ay, bx, by);
                        /* Compute squared distance between points */
                        /* 'C' and 'D'.                            */
                        double cd2 = SquaredDistance(cx, cy, dx, dy);
                        /* If CD distance is not near zero */
                        /* (ex. 0.5) ...                   */
                        if (cd2 > lqmParams.pores_min_dist2) {
                            /* Compute ratio of squared distances. */
                            double ratio = ab2 / cd2;

                            /* If ratio is small enough (ex. 2.25)...*/
                            if (ratio <= lqmParams.pores_max_ratio) {
                                /* Then assume pore & remove minutia. */
                                minutiae.erase(minutiae.begin() + i);
                                continue;
                            }
                            /* Otherwise, ratio to big, so assume */
                            /* legitimate minutia.                */
                        } /* Else, cd2 too small. */
                    } /* Else, R is on pixel the same color as type, so do not */
                    /* remove minutia point and skip to next one.            */
                } /* Else block is unreliable or has INVALID direction. */

                /* Current minutia was not removed. */
                /* Bump to next minutia in list. */
                ++i;
                /* Otherwise, next minutia has slid into slot of current removed one. */
            } /* End While minutia remaining in list. */
        }

        bool SearchInDirection(FeaturePoint& op, const int pix,
                        const int strt_x, const int strt_y,
                        const double delta_x, const double delta_y, const int maxsteps,
                        unsigned char *bdata, const int iw, const int ih)
        {
            /* Set previous point to starting point. */
            int px = strt_x;
            int py = strt_y;
            /* Set floating point accumulators to starting point. */
            double fx = static_cast<double>(strt_x);
            double fy = static_cast<double>(strt_y);

            /* Foreach step up to the specified maximum ... */
            for (int i = 0; i < maxsteps; ++i) {

                /* Increment accumulators. */
                fx += delta_x;
                fy += delta_y;
                /* Round to get next step. */
                int x = sround(fx);
                int y = sround(fy);

                /* If we stepped outside the image boundaries ... */
                if (x < 0 || x >= iw || y < 0 || y >= ih) {
                    break;
                }

                /* Otherwise, test to see if we found our pixel with value 'pix'. */
                if (*(bdata+(y*iw)+x) == pix) {
                    /* The previous and current pixels form a feature, edge pixel */
                    /* pair, which we would like to use for edge following.  The  */
                    /* previous pixel may be a diagonal neighbor however to the   */
                    /* current pixel, in which case the pair could not be used by */
                    /* the contour tracing (which requires the edge pixel in the  */
                    /* pair neighbor to the N,S,E or W.                           */
                    /* This routine adjusts the pair so that the results may be   */
                    /* used by the contour tracing.                               */
                    op.x = x;
                    op.y = y;
                    op.ex = px;
                    op.ey = py;
                    FixEdgePixelPair(op, bdata, iw);

                    /* Return true (we found what we were looking for). */
                    return true;
                }

                /* Otherwise, still haven't found pixel with desired value, */
                /* so set current point to previous and take another step.  */
                px = x;
                py = y;
            }

            /* Return false (we did not find what we were looking for). */
            op.x = -1;
            op.y = -1;
            op.ex = -1;
            op.ey = -1;
            return false;
        }

        void FixEdgePixelPair(FeaturePoint& p, unsigned char *bdata, const int iw) {
            /* Get the pixel value of the feature. */
            int feature_pix = *(bdata + (p.y * iw) + p.x);

            /* Store the input points to current and previous points. */
            int cx = p.x;
            int cy = p.y;
            int px = p.ex;
            int py = p.ey;

            /* Compute detlas between current and previous point. */
            int dx = px - cx;
            int dy = py - cy;

            /* If previous point (P) is diagonal neighbor of    */
            /* current point (C)... This is a problem because   */
            /* the contour tracing routine requires that the    */
            /* "edge" pixel be north, south, east, or west of   */
            /* of the feature point.  If the previous pixel is  */
            /* diagonal neighbor, then we need to adjust either */
            /* the positon of te previous or current pixel.     */
            if (std::abs(dx) == 1 && std::abs(dy) == 1) {
                /* Then we have one of the 4 following conditions:  */
                /*                                                  */
                /*        *C                             C*         */
                /*    1.  P*     2.  P*     3. *P     4. *P         */
                /*                   *C        C*                   */
                /*                                                  */
                /*  dx =  -1         -1         1         1         */
                /*  dy =   1         -1        -1         1         */
                /*                                                  */
                /* Want to test values in positions of '*':         */
                /*  Let point P == (px, py)                         */
                /*           p1 == '*' positon where x changes      */
                /*           p2 == '*' positon where y changes      */
                /*                                                  */
                /*  p1 = px+1,py    px+1,py   px-1,py    px-1,py    */
                /*  p2 = px,py-1    px,py+1   px,py+1    px,py-1    */
                /*                                                  */
                /* These can all be rewritten:                      */
                /*  p1 = px-dx,py                                   */
                /*  p2 = px,py-dy                                   */

                /* Check if 'p1' is NOT the value we are searching for... */
                if (*(bdata+(py*iw)+(px-dx)) != feature_pix) {
                    /* Then set x-coord of edge pixel to p1. */
                    px -= dx;
                }
                /* Check if 'p2' is NOT the value we are searching for... */
                else if (*(bdata+((py-dy)*iw)+px) != feature_pix) {
                    /* Then set y-coord of edge pixel to p2. */
                    py -= dy;
                }
                /* Otherwise, the current pixel 'C' is exposed on a corner ... */
                else {
                    /* Set pixel 'C' to 'p1', which also has the pixel */
                    /* value we are searching for.                     */
                    cy += dy;
                }

                /* Set the pointers to the resulting values. */
                p.x = cx;
                p.y = cy;
                p.ex = px;
                p.ey = py;
            }

            /* Otherwise, nothing has changed. */
        }

        void CountMinutiaeRidges(Minutiae& minutiae,
            unsigned char *bdata, const int iw, const int ih,
            const LQMParams& lqmParams)
        {
            /* Sort minutia points on x then y (column-oriented). */
            SortMinutiaeXY(minutiae, iw);

            /* Remove any duplicate minutia points from the list. */
            RemoveDuplicateMinutiae(minutiae);

            /* Foreach remaining sorted minutia in list ... */
            for (int i = 0; i < (static_cast<int>(minutiae.size())) - 1; ++i) {
                /* Located neighbors and count number of ridges in between. */
                /* NOTE: neighbor and ridge count results are stored in     */
                /*       minutiae[i]                                        */
                CountMinutiaRidges(i, minutiae, bdata, iw, ih, lqmParams);
            }
        }

        void CountMinutiaRidges(const int first, Minutiae& minutiae,
                            unsigned char *bdata, const int iw, const int ih,
                            const LQMParams& lqmParams)
        {
            /* Find up to the maximum number of qualifying neighbors. */
            std::vector<int> nbr_list;
            FindNeighbors(nbr_list, lqmParams.max_nbrs, first, minutiae);

            /* If no neighors found ... */
            if (nbr_list.size() == 0) {
                /* Then no list returned and no ridges to count. */
                return;
            }

            /* Sort neighbors on delta dirs. */
            SortNeighbors(nbr_list, first, minutiae);

            /* Count ridges between first and neighbors. */
            /* List of ridge counts, one for each neighbor stored. */
            std::vector<int> nbr_nridges;
            nbr_nridges.reserve(nbr_list.size());

            /* Foreach neighbor found and sorted in list ... */
            for (int i = 0; i < static_cast<int>(nbr_list.size()); ++i) {
                /* Count the ridges between the primary minutia and the neighbor. */
                int count = RidgeCount(first, nbr_list[static_cast<std::size_t>(i)], minutiae, bdata, iw, ih, lqmParams);

                /* Otherwise, ridge count successful, so store ridge count to list. */
                nbr_nridges.push_back(count);
            }

            /* Assign neighbor indices and ridge counts to primary minutia. */
            minutiae[static_cast<std::size_t>(first)].nbrs = nbr_list;
            minutiae[static_cast<std::size_t>(first)].ridge_counts = nbr_nridges;
        }

        void FindNeighbors(std::vector<int>& onbr_list, const int max_nbrs, const int first, Minutiae& minutiae)
        {
            /* Allocate list of neighbor minutiae indices. */
            onbr_list.reserve(static_cast<std::size_t>(max_nbrs));

            /* Allocate list of squared euclidean distances between neighbors */
            /* and current primary minutia point.                             */
            std::vector<double> nbr_sqr_dists;
            nbr_sqr_dists.reserve(static_cast<std::size_t>(max_nbrs));

            /* Assign secondary to one passed current primary minutia. */
            int second = first + 1;
            /* Compute location of maximum last stored neighbor. */
            int last_nbr = max_nbrs - 1;

            /* While minutia (in sorted order) still remian for processing ... */
            /* NOTE: The minutia in the input list have been sorted on X and   */
            /* then on Y.  So, the neighbors are selected according to those   */
            /* that lie below the primary minutia in the same pixel column and */
            /* then subsequently those that lie in complete pixel columns to   */
            /* the right of the primary minutia.                               */
            while (second < static_cast<int>(minutiae.size())) {
                /* Assign temporary minutia pointers. */
                const Minutia& minutia1 = minutiae[static_cast<std::size_t>(first)];
                const Minutia& minutia2 = minutiae[static_cast<std::size_t>(second)];

                /* Compute squared distance between minutiae along x-axis. */
                double xdist = minutia2.x - minutia1.x;
                double xdist2 = xdist * xdist;

                /* If the neighbor lists are not full OR the x-distance to current */
                /* secondary is smaller than maximum neighbor distance stored ...  */
                if (static_cast<int>(onbr_list.size()) < max_nbrs || xdist2 < nbr_sqr_dists[static_cast<std::size_t>(last_nbr)]) {
                    /* Append or insert the new neighbor into the neighbor lists. */
                    UpdateNeighborDists(onbr_list, nbr_sqr_dists, max_nbrs, first, second, minutiae);
                }
                /* Otherwise, if the neighbor lists is full AND the x-distance   */
                /* to current secondary is larger than maximum neighbor distance */
                /* stored ...                                                    */
                else {
                    /* So, stop searching for more neighbors. */
                    break;
                }

                /* Bump to next secondary minutia. */
                ++second;
            }
        }

        void UpdateNeighborDists(std::vector<int>& nbr_list, std::vector<double>& nbr_sqr_dists,
            const int max_nbrs, const int first, const int second, Minutiae& minutiae)
        {
            /* Compute position of maximum last neighbor stored. */
            int last_nbr = max_nbrs - 1;

            /* Assigne temporary minutia pointers. */
            const Minutia& minutia1 = minutiae[static_cast<std::size_t>(first)];
            const Minutia& minutia2 = minutiae[static_cast<std::size_t>(second)];

            /* Compute squared euclidean distance between minutia pair. */
            double dist2 = SquaredDistance(minutia1.x, minutia1.y, minutia2.x, minutia2.y);

            /* If maximum number of neighbors not yet stored in lists OR */
            /* if the squared distance to current secondary is less      */
            /* than the largest stored neighbor distance ...             */
            if (static_cast<int>(nbr_list.size()) < max_nbrs || dist2 < nbr_sqr_dists[static_cast<std::size_t>(last_nbr)]) {
                /* Find insertion point in neighbor lists. */
                int pos = FindIncreasingPositionDouble(dist2, nbr_sqr_dists);
                /* If the position returned is >= maximum list length (this should */
                /* never happen, but just in case) ...                             */
                if (pos >= max_nbrs) {
                    throw std::invalid_argument("The position index returned from FindIncreasingPositionDouble for new neighbor exceeds max_nbrs");
                }
                /* Insert the new neighbor into the neighbor lists at the */
                /* specified location.                                    */
                InsertNeighbor(pos, second, dist2, nbr_list, nbr_sqr_dists, max_nbrs);

                /* Otherwise, neighbor inserted successfully, so return */
                return;
            }
            /* Otherwise, the new neighbor is not sufficiently close to be       */
            /* added or inserted into the neighbor lists, so ignore the neighbor */
        }

        int FindIncreasingPositionDouble(const double val, std::vector<double>& list) {
            int i = 0;
            /* Foreach item in double list ... */
            for (; i < static_cast<int>(list.size()); ++i) {
                /* If the value is smaller than the current item in list ... */
                if (val < list[static_cast<std::size_t>(i)]) {
                    /* Then we found were to insert the value in the list maintaining */
                    /* an increasing sorted order.                                    */
                    return i;
                }

                /* Otherwise, the value is still larger than current item, so */
                /* continue to next item in the list.                         */
            }

            /* Otherwise, we never found a slot within the list to insert the */
            /* the value, so place at the end of the sorted list.             */
            return i;
        }

        void InsertNeighbor(const int pos, const int nbr_index, const double nbr_dist2,
            std::vector<int>& nbr_list, std::vector<double>& nbr_sqr_dists, const int max_nbrs)
        {
            /* If the desired insertion position is beyond one passed the last     */
            /* neighbor in the lists OR greater than equal to the maximum ...      */
            /* NOTE: pos is zero-oriented while nbr_list.size() and max_nbrs are 1-oriented. */
            if (pos > static_cast<int>(nbr_list.size()) || pos >= max_nbrs) {
                throw std::invalid_argument("InsertNeighbor invalid argument: insertion point is out of bounds");
            }

            /* If the neighbors lists are full ... */
            if (static_cast<int>(nbr_list.size()) == max_nbrs) {
                /* So, we must bump the last neighbor in the lists off to make */
                /* room for the new neighbor (ignore last neighbor in lists).  */
                nbr_list.pop_back();
                nbr_sqr_dists.pop_back();
            }
            /* Else if there is a list overflow error condition */
            /* (shouldn't ever happen, but just in case) ...       */
            else if (static_cast<int>(nbr_list.size()) > max_nbrs) {
                throw std::invalid_argument("InsertNeighbor invalid argument: size of nbr_list exceeds max_nbrs");
            }

            /* Insert the new neighbor into the lists at pos */
            nbr_list.insert(nbr_list.begin() + pos, nbr_index);
            nbr_sqr_dists.insert(nbr_sqr_dists.begin() + pos, nbr_dist2);
        }

        void SortNeighbors(std::vector<int>& nbr_list, const int first, Minutiae& minutiae) {
            /* List of angles of lines joining the current primary to each */
            /* of the secondary neighbors.                                 */
            std::vector<std::pair<double, int>> join_thetas;
            join_thetas.reserve(nbr_list.size());

            for (int i = 0; i < static_cast<int>(nbr_list.size()); ++i) {
            /* Compute angle to line connecting the 2 points.             */
            /* Coordinates are swapped and order of points reversed to    */
            /* account for 0 direction is vertical and positive direction */
            /* is clockwise.                                              */
            double theta = AngleToLine(minutiae[static_cast<std::size_t>(nbr_list[static_cast<std::size_t>(i)])].y,
                                minutiae[static_cast<std::size_t>(nbr_list[static_cast<std::size_t>(i)])].x,
                                minutiae[static_cast<std::size_t>(first)].y,
                                minutiae[static_cast<std::size_t>(first)].x);

            /* Make sure the angle is positive. */
            theta += LQM_PI2;
            theta = std::fmod(theta, LQM_PI2);
            join_thetas.emplace_back(theta, nbr_list[static_cast<std::size_t>(i)]);
            }

            /* Sort the neighbor indicies into rank order. */
            std::sort(join_thetas.begin(), join_thetas.end(), [&](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                return a.first < b.first;
            });

            /* Copy sorted indices from join_thetas to nbr_list */
            for (int i = 0; i < static_cast<int>(join_thetas.size()); ++i) {
                nbr_list[static_cast<std::size_t>(i)] = join_thetas[static_cast<std::size_t>(i)].second;
            }
        }

        int RidgeCount(const int first, const int second, Minutiae& minutiae,
            unsigned char *bdata, const int iw, const int ih,
            const LQMParams& lqmParams)
        {
            const Minutia& minutia1 = minutiae[static_cast<std::size_t>(first)];
            const Minutia& minutia2 = minutiae[static_cast<std::size_t>(second)];

            /* If the 2 mintuia have identical pixel coords ... */
            if (minutia1.x == minutia2.x && minutia1.y == minutia2.y) {
                /* Then zero ridges between points. */
                return 0;
            }

            /* Compute linear trajectory of contiguous pixels between first */
            /* and second minutia points.                                   */
            std::vector<std::pair<int, int>> points;
            OpenLQM::Core::LinePoints(points, minutia1.x, minutia1.y, minutia2.x, minutia2.y);

            /* It there are no points on the line trajectory, then no ridges */
            /* to count (this should not happen, but just in case) ...       */
            if (points.size() == 0) {
                return 0;
            }

            /* Find first pixel opposite type along linear trajectory from */
            /* first minutia.                                              */
            int prevpix = *(bdata + (points[0].second*iw) + points[0].first);
            int i = 1;
            bool found = false;
            while (i < static_cast<int>(points.size())) {
            int curpix = *(bdata + (points[static_cast<std::size_t>(i)].second*iw) + points[static_cast<std::size_t>(i)].first);
            if (curpix != prevpix) {
                found = true;
                break;
            }
            ++i;
            }

            /* If opposite pixel not found ... then no ridges to count */
            if (!found) {
                return 0;
            }

            /* Ready to count ridges, so initialize counter to 0. */
            int ridge_count = 0;

            /* While not at the end of the trajectory ... */
            while (i < static_cast<int>(points.size())) {
                /* If 0-to-1 transition not found ... */
                if (!FindTransition(&i, 0, 1, points, bdata, iw)) {
                        /* Then we are done looking for ridges. */
                        /* Return number of ridges counted to this point. */
                        return ridge_count;
                }
                /* Otherwise, we found a new ridge start transition, so store */
                /* its location (the location of the 1 in 0-to-1 transition). */
                int ridge_start = i;

                /* If 1-to-0 transition not found ... */
                if (!FindTransition(&i, 1, 0, points, bdata, iw)) {
                    /* Then we are done looking for ridges. */
                    /* Return number of ridges counted to this point. */
                    return ridge_count;
                }
                /* Otherwise, we found a new ridge end transition, so store   */
                /* its location (the location of the 0 in 1-to-0 transition). */
                int ridge_end = i;

                /* Conduct the validation, tracing the contour of the ridge  */
                /* from the ridge ending point a specified number of steps   */
                /* scanning for neighbors clockwise and counter-clockwise.   */
                /* If the ridge starting point is encounted during the trace */
                /* then we can assume we do not have a valid ridge crossing  */
                /* and instead we are walking on and off the edge of the     */
                /* side of a ridge.                                          */
                bool validated = ValidateRidgeCrossing(ridge_start, ridge_end,
                                                points, bdata, iw, ih,
                                                lqmParams.max_ridge_steps);

                /* If validation result is TRUE ... */
                if (validated) {
                    /* Then assume we have found a valid ridge crossing and bump */
                    /* the ridge counter.                                        */
                    ++ridge_count;
                }

                /* Otherwise, ignore the current ridge start and end transitions */
                /* and go back and search for new ridge start.                   */
            }

            /* Return the number of ridges counted. */
            return ridge_count;
        }

        bool FindTransition(int *iptr, const int pix1, const int pix2,
            const std::vector<std::pair<int, int>>& points,
            unsigned char *bdata, const int iw)
        {
            /* Set previous index to starting position. */
            int i = *iptr;
            /* Bump previous index by 1 to get next index. */
            int j = i+1;

            /* While not one point from the end of the trajectory .. */
            while (i < (static_cast<int>(points.size())) - 1) {
                /* If we have found the desired transition ... */
                if ((*(bdata+(points[static_cast<std::size_t>(i)].second*iw) + points[static_cast<std::size_t>(i)].first) == pix1) && (*(bdata+(points[static_cast<std::size_t>(j)].second*iw) + points[static_cast<std::size_t>(j)].first) == pix2)) {
                    /* Adjust the position pointer to the location of the */
                    /* second pixel in the transition.                    */
                    *iptr = j;

                    /* Return TRUE. */
                    return true;
                }
                /* Otherwise, the desired transition was not found in current */
                /* pixel pair, so bump to the next pair along the trajector.  */
                ++i;
                ++j;
            }

            /* If we get here, then we exhausted the trajector without finding */
            /* the desired transition, so set the position pointer to the end  */
            /* of the trajector, and return FALSE.                             */
            *iptr = static_cast<int>(points.size());
            return false;
        }

        bool ValidateRidgeCrossing(const int ridge_start, const int ridge_end,
            const std::vector<std::pair<int, int>>& points,
            unsigned char *bdata, const int iw, const int ih,
            const int max_ridge_steps)
        {
            /* Assign edge pixel pair for contour trace. */
            FeaturePoint feat;
            feat.x = points[static_cast<std::size_t>(ridge_end)].first;
            feat.y = points[static_cast<std::size_t>(ridge_end)].second;
            feat.ex = points[static_cast<std::size_t>(ridge_end-1)].first;
            feat.ey = points[static_cast<std::size_t>(ridge_end-1)].second;

            /* Adjust pixel pair if they neighbor each other diagonally. */
            FixEdgePixelPair(feat, bdata, iw);

            /* Trace ridge contour, starting at the ridge end transition, and */
            /* taking a specified number of step scanning for edge neighbors  */
            /* clockwise.  As we trace the ridge, we want to detect if we     */
            /* encounter the ridge start transition.  NOTE: The ridge end     */
            /* position is on the white (of a black to white transition) and  */
            /* the ridge start is on the black (of a black to white trans),   */
            /* so the edge trace needs to look for the what pixel (not the    */
            /* black one) of the ridge start transition.                      */
            int ret;
            {
                Contour contour;
                ret = TraceContour(contour,
                                    max_ridge_steps,
                                    points[static_cast<std::size_t>(ridge_start-1)].first, points[static_cast<std::size_t>(ridge_start-1)].second,
                                    feat.x, feat.y, feat.ex, feat.ey,
                                    LQM_SCAN_CLOCKWISE, bdata, iw, ih);
            }

            if (ret == LQM_IGNORE) {
                /* Some sort of initialization error, return false. */
                return false;
            }

            /* Treat this the same as if the trace were actually located at the  */
            /* ridge start point (in which case LOOP_FOUND is returned).         */
            /* So, If not IGNORED and ridge start not encounted in trace ...     */
            if (ret != LQM_LOOP_FOUND) {
                /* Now conduct contour trace scanning for edge neighbors counter- */
                /* clockwise.                                                     */
                Contour contour;
                ret = TraceContour(contour,
                                    max_ridge_steps,
                                    points[static_cast<std::size_t>(ridge_start-1)].first, points[static_cast<std::size_t>(ridge_start-1)].second,
                                    feat.x, feat.y, feat.ex, feat.ey,
                                    LQM_SCAN_COUNTER_CLOCKWISE, bdata, iw, ih);

                /* Otherwise, if the trace was not IGNORED, then a contour was */
                /* was generated and returned.                                 */

                /* If trace not IGNORED and ridge start not encounted in 2nd trace ... */
                if (ret != LQM_IGNORE && ret != LQM_LOOP_FOUND) {
                    /* If we get here, assume we have a ridge crossing. */
                    return true;
                }
                /* Otherwise, second trace returned IGNORE or ridge start found. */
            }
            /* Otherwise, first trace returned IGNORE or ridge start found. */

            /* If we get here, then we failed to validate a ridge crossing. */
            return false;
        }

        void RemoveDuplicateMinutiae(Minutiae& minutiae) {
            /* Work backward from the end of the list of minutiae.  This way */
            /* we can selectively remove minutia from the list and not cause */
            /* problems with keeping track of current indices.               */
            int i = (static_cast<int>(minutiae.size())) - 1;

            while (i > 0) {
                const Minutia& minutia1 = minutiae[static_cast<std::size_t>(i)];
                const Minutia& minutia2 = minutiae[static_cast<std::size_t>(i - 1)];
                /* If minutia pair has identical coordinates ... */
                if (minutia1.x == minutia2.x && minutia1.y == minutia2.y) {
                    /* Remove the 1st minutia from the minutiae list. */
                    minutiae.erase(minutiae.begin() + i - 1);
                }

                /* Decrement to compare next minutia pair in list. */
                --i;
            }
        }

        int FilterAndOutputMinutiae(MinutiaOut* pOut, int maxMinutiae, Minutiae& minutiae, int inputResolution) {
            const int offset = 3;
            int resolutionFactor = inputResolution / 500;

            int maxIndex = std::min(maxMinutiae, static_cast<int>(minutiae.size()));
            MinutiaOut* pOutMinutia = pOut;
            for (int i = 0; i < maxIndex; ++i) {
                Minutia& minutia = minutiae[static_cast<std::size_t>(i)];

                if (FindCloseMinutiae(minutiae, minutia, 4)) {
                    // Not updating maxIndex in order to match original behavior of truncating minutiae beyond maxMinutiae prior to
                    //     filtering with FindCloseMinutiae
                    continue;
                }

                pOutMinutia->x = minutia.x * resolutionFactor;
                pOutMinutia->y = minutia.y * resolutionFactor;
                pOutMinutia->type = minutia.type;
                pOutMinutia->theta = AngleToDegrees(minutia.direction);

                int absoluteOffset = -offset * resolutionFactor;

                std::pair<int, int> correctedPoint = GetDirectedPoint(pOutMinutia->x, pOutMinutia->y, pOutMinutia->theta, absoluteOffset);
                pOutMinutia->x = correctedPoint.first;
                pOutMinutia->y = correctedPoint.second;

                ++pOutMinutia;
            }

            return static_cast<int>(pOutMinutia - pOut);
        }

        bool FindCloseMinutiae(const Minutiae& minutiae, const Minutia& testMinutia, int pixelDistance) {
        /*
            - Assumes minutiae are sorted by X
            - Assumes features are never placed within +-pixelDistance pixels of each other
            - pixelDistance is always in 500dpi pixels
        */

            for (const Minutia& minutia : minutiae) {
                if (&minutia == &testMinutia) {
                    continue;
                }

                if (minutia.x > (testMinutia.x - pixelDistance)) {
                    if (minutia.x > (testMinutia.x + pixelDistance)) {
                        // We've moved past all possible candidates
                        return false;
                    } else if (minutia.y > (testMinutia.y - pixelDistance) && minutia.y < (testMinutia.y + pixelDistance)) {
                        return true;
                    }
                }
            }

            return false;
        }

        long AngleToDegrees(int angle) {
            /*
                Directions 1..32 each rotates CCW 1/32 of a circle.
                When angle is 0, the direction is 90 degrees. Therefore, we have to
                    divide 360 degrees into 32 steps, then multiply by the angle and add 90 degrees.
                The incoming angle also needs to be inverted, hence subtracting it from 360.
            */

            // Casting double to float here is done to maintain conformance with the original output
            float ret = static_cast<float>((360.0 - (static_cast<double>(angle) * (360.0/32.0))) + 90.0);
            while (ret < 0.0f) {
                ret += 360.0f;
            }

            while (ret >= 360.0f) {
                ret -= 360.0f;
            }

            // This rounding function is used to match the behavior of Visual Basic's
            //  conversion from floating point types to integral types (banker's rounding)
            return std::lrint(ret);
        }

        std::pair<int, int> GetDirectedPoint(int x, int y, double angle, double length) {
            // This double -> float narrowing conversion was ported as-is to maintain output conformance
            float angleRads = static_cast<float>(angle * (LQM_VB_PI / 180.0));

            // The float -> double widening conversions of the trig function arguments were ported as-is to maintain output conformance
            return std::pair<int, int>(x + std::lrint(length * cos(static_cast<double>(angleRads))), y + std::lrint(length * sin(static_cast<double>(angleRads))));
        }

        void FilterMinutiaeByQuality(const unsigned char* pLocQ, int mapWidth, int mapHeight, int inputWidth, int inputHeight, MinutiaOut* pMinutiae, int minutiaCount, std::vector<MinutiaOut>& sortedMinutiae) {
            for (MinutiaOut* pMinutia = pMinutiae; pMinutia < pMinutiae + minutiaCount; ++pMinutia) {
                int xQuality = static_cast<int>(std::lrint(static_cast<double>(pMinutia->x) * (static_cast<double>(mapWidth) / static_cast<double>(inputWidth))));
                int yQuality = static_cast<int>(std::lrint(static_cast<double>(pMinutia->y) * (static_cast<double>(mapHeight) / static_cast<double>(inputHeight))));

                unsigned char quality;

                if (xQuality < 0 || yQuality < 0) {
                    quality = 0;
                } else {
                    xQuality = std::min(xQuality, mapWidth - 1);
                    yQuality = std::min(yQuality, mapHeight - 1);
                    quality = pLocQ[yQuality * mapWidth + xQuality];
                }

                if (quality >= 2) {
                    sortedMinutiae.emplace_back(pMinutia->x, pMinutia->y, pMinutia->type, static_cast<long>(quality));
                }
            }

            std::sort(sortedMinutiae.begin(), sortedMinutiae.end(), [](const MinutiaOut& a, const MinutiaOut& b) {
                return a.theta < b.theta;
            });
        }

        void CalculateMinutiaCounts(const std::vector<MinutiaOut>& minutiae, std::array<int, 3>& minutiaCounts) {
            for (int& count : minutiaCounts) {
                count = 0;
            }

            for (const MinutiaOut& minutia : minutiae) {
                int qualityIndex = minutia.theta - 2;
                if (qualityIndex >= 0 && qualityIndex < static_cast<int>(minutiaCounts.size())) {
                    ++minutiaCounts[static_cast<std::size_t>(qualityIndex)];
                }
            }
        }
    }
}
