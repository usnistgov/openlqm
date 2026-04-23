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

#pragma once

#include <vector>
#include <utility>

// todo: move these mathematical definitions elsewhere
#define LQM_PI 3.1415926535897932
#define LQM_PI2 (LQM_PI * 2.0)

// This alternative value of PI is used in code ported from Visual Basic
#define LQM_VB_PI 3.1415926535897931

#define sround(x) (static_cast<int> (((x)<0) ? (x)-0.5 : (x)+0.5))
/* Truncate floating point precision by multiply, rounding, and then */
/* dividing by this value.  This enables consistant results across   */
/* different computer architectures.                                 */
#define TRUNC_SCALE          16384.0
#define trunc_dbl_precision(x) (static_cast<double>(((x)<0.0) \
                 ? (static_cast<int>(((x)*(TRUNC_SCALE))-0.5))/(TRUNC_SCALE) \
                 : (static_cast<int>(((x)*(TRUNC_SCALE))+0.5))/(TRUNC_SCALE)))




/***** DFT CONSTANTS *****/

/* This specifies the number of DFT wave forms to be applied */
#define NUM_DFT_WAVES            4

/* Various thresholds and factors used by OpenLQM.   */

/* Minimum DFT power allowable in any one direction. */
#define POWMAX_MIN          100000.0

/* Minimum normalized power allowable in any one     */
/* direction.                                        */
#define POWNORM_MIN              3.8

/* Maximum power allowable at the lowest frequency   */
/* DFT wave.                                         */
#define POWMAX_MAX        50000000.0

/* Check for a fork at +- this number of units from  */
/* current integer direction.  For example,          */
/*           2 dir ==> 11.25 X 2 degrees.            */
#define FORK_INTERVAL            2  * NUM_DIRECTIONS_FACTOR

/* Minimum DFT power allowable at fork angles is     */
/* FORK_PCT_POWMAX X block's max directional power.  */
#define FORK_PCT_POWMAX          0.7

/* Minimum normalized power allowable at fork angles */
/* is FORK_PCT_POWNORM X POWNORM_MIN                 */
#define FORK_PCT_POWNORM         0.75

/* Designates LQM parameter as undefined. */
#define LQM_UNDEFINED               -1

/* Minimum strength for a direction to be considered significant. */
#define LQM_DIR_STRENGTH_MIN         0.2

/* Definitions for 8-bit binary pixel intensities. */
#define LQM_WHITE_PIXEL            255
#define LQM_BLACK_PIXEL              0


/* Initial reserve capacity of the list of detected minutia  */
#define LQM_INITIAL_MINUTIA_CAPACITY          1000

/*************************************************************************/
/* 10, 2X3 pixel pair feature patterns used to define ridge endings      */
/* and bifurcations.                                                     */
/* 2nd pixel pair is permitted to repeat multiple times in match.        */
#define LQM_NFEATURES      10
#define LQM_BIFURCATION     0
#define LQM_RIDGE_ENDING    1
#define LQM_DISAPPEARING    0
#define LQM_APPEARING       1

// In EmbedMinutiaeDetail, look for neighbors this many blocks away:
#define LQM_NEIGHBOR_DELTA 2

/* Minutia scanning constants */

// Returned from OnHook if a pair of minutiae lie on a hook on the side of a ridge or valley
#define LQM_HOOK_FOUND 1
// Returned from minutia scanning functions when a contour loop is found
#define LQM_LOOP_FOUND 1
// Returned from minutia scanning functions when minutia should be ignored
#define LQM_IGNORE 2
// Returned from GetCenteredContour when either half-contour is incomplete
#define LQM_INCOMPLETE 3

/* Map value designating a block is near a high-curvature */
/* area such as a core or delta.                          */
#define LQM_HIGH_CURVATURE -2

/* Definitions for controlling the scanning of minutiae. */
#define LQM_SCAN_HORIZONTAL          0
#define LQM_SCAN_VERTICAL            1
#define LQM_SCAN_CLOCKWISE           0
#define LQM_SCAN_COUNTER_CLOCKWISE   1

namespace OpenLQM {
   namespace Core {
      enum GridOffsetMode {
         /* Designates that rotated grid offsets should be relative */
         /* to the grid's center.                                   */
         RELATIVE2CENTER = 0,
         /* Designates that rotated grid offsets should be relative */
         /* to the grid's origin.                                   */
         RELATIVE2ORIGIN = 1
      };

      struct FeaturePattern {
         int type;
         int appearing;
         int first[2];
         int second[2];
         int third[2];
      };

      /* Parameters used by LQM for setting thresholds and  */
      /* defining testing criterion.                        */
      struct LQMParams {
         /* Image Controls */
         int    pad_value;
         int    join_line_radius;

         /* Map Controls */
         int    blocksize;       /* Pixel dimension image block.                 */
         int    windowsize;      /* Pixel dimension window surrounding block.    */
         int    windowoffset;    /* Offset in X & Y from block to window origin. */
         int    num_directions;
         double start_dir_angle;
         int    rmv_valid_nbr_min;
         double dir_strength_min;
         int    dir_distance_max;
         int    smth_valid_nbr_min;
         int    vort_valid_nbr_min;
         int    highcurv_vorticity_min;
         int    highcurv_curvature_min;
         int    min_interpolate_nbrs;
         int    percentile_min_max;
         int    min_contrast_delta;

         /* DFT Controls */
         int    num_dft_waves;
         double powmax_min;
         double pownorm_min;
         double powmax_max;
         int    fork_interval;
         double fork_pct_powmax;
         double fork_pct_pownorm;

         /* Binarization Controls */
         int    dirbin_grid_w;
         int    dirbin_grid_h;
         int    isobin_grid_dim;
         int    num_fill_holes;

         /* Minutiae Detection Controls */
         int    max_minutia_delta;
         double max_high_curve_theta;
         int    high_curve_half_contour;
         int    min_loop_len;
         double min_loop_aspect_dist;
         double min_loop_aspect_ratio;

         /* Minutiae Link Controls */
         int    link_table_dim;
         int    max_link_dist;
         int    min_theta_dist;
         int    maxtrans;
         double score_theta_norm;
         double score_dist_norm;
         double score_dist_weight;
         double score_numerator;

         /* False Minutiae Removal Controls */
         int    max_rmtest_dist;
         int    max_hook_len;
         int    max_half_loop;
         int    trans_dir_pix;
         int    small_loop_len;
         int    side_half_contour;
         int    inv_block_margin;
         int    rm_valid_nbr_min;
         int    max_overlap_dist;
         int    max_overlap_join_dist;
         int    malformation_steps_1;
         int    malformation_steps_2;
         double min_malformation_ratio;
         int    max_malformation_dist;
         int    pores_trans_r;
         int    pores_perp_steps;
         int    pores_steps_fwd;
         int    pores_steps_bwd;
         double pores_min_dist2;
         double pores_max_ratio;

         /* Ridge Counting Controls */
         int    max_nbrs;
         int    max_ridge_steps;
      };

      struct Minutia {
         int x;
         int y;
         int ex;
         int ey;
         int direction;
         double reliability;
         int type;
         int appearing;
         int feature_id;
         std::vector<int> nbrs;
         std::vector<int> ridge_counts;

         Minutia(
            const int x_loc, const int y_loc,
            const int x_edge, const int y_edge, const int idir,
            const double reliability_,
            const int type_, const int appearing_, const int feature_id_
         ) {
            /* Assign minutia structure attributes. */
            this->x = x_loc;
            this->y = y_loc;
            this->ex = x_edge;
            this->ey = y_edge;
            this->direction = idir;
            this->reliability = reliability_;
            this->type = type_;
            this->appearing = appearing_;
            this->feature_id = feature_id_;
         }
      };

      typedef std::vector<Minutia> Minutiae;

      struct Row {
         int y;                  /* Y-coord of current row in shape.                  */
         std::vector<int> xs;    /* X-coords for shape contour points on current row. */

         Row(int y_, int alloc_pts) {
            this->y = y_;
            this->xs.reserve(static_cast<std::size_t>(alloc_pts));
         }
      };

      struct Shape {
         int ymin;      /* Y-coord of top-most scanline in shape.     */
         int ymax;      /* Y-coord of bottom-most scanline in shape.  */
         std::vector<Row> rows; /* List of rows comprising the shape. */

         Shape(const int *contour_x, const int *contour_y, const int ncontour);
      };

      enum ExtremumType {
         MINIMUM = -1,
         INDETERMINATE = 0,
         MAXIMUM = 1
      };

      struct Extremum {
         int val;
         ExtremumType type;
         int sourceIndex;

         Extremum(int val_, ExtremumType type_, int sourceIndex_) {
            this->val = val_;
            this->type = type_;
            this->sourceIndex = sourceIndex_;
         }
      };

      struct MinutiaOut {
         long x;
         long y;
         long type;
         long theta;

         MinutiaOut() {

         }

         MinutiaOut(long x_, long y_, long type_, long theta_) {
            x = x_;
            y = y_;
            type = type_;
            theta = theta_;
         }
      };

      /* Lookup tables for converting from integer directions */
      /* to angles in radians.                                */
      struct Dir2Rad {
         int ndirs;
         std::vector<double> cosLut;
         std::vector<double> sinLut;

         explicit Dir2Rad(int ndirs_);
      };

      /* DFT wave form structure containing both cosine and   */
      /* sine components for a specific frequency.            */
      struct DFTWave {
         std::vector<double> cosLut;
         std::vector<double> sinLut;

         void Init(int blocksize, double freq);
      };

      /* DFT wave forms structure containing all wave forms  */
      /* to be used in DFT analysis.                         */
      struct DFTWaves {
         int nwaves;
         int wavelen;
         std::vector<DFTWave> waves;

         DFTWaves(int nwaves_, int blocksize);
      };

      /* Rotated pixel offsets for a grid of specified dimensions */
      /* rotated at a specified number of different orientations  */
      /* (directions).  This structure used by the DFT analysis   */
      /* when generating a Direction Map and also for conducting  */
      /* isotropic binarization.                                  */
      struct RotGrids {
         int pad;
         GridOffsetMode relative2;
         double start_angle;
         int ngrids;
         int grid_w;
         int grid_h;
         std::vector<std::vector<int>> grids;

      /*************************************************************************
       Constructor     - Allocates and initializes a set of offsets that address
                     individual rotated pixels within a grid.
                     These rotated grids are used to conduct DFT analyses
                     on blocks of input image data, and they are used
                     in isotropic binarization.

         Input:
            iw        - width (in pixels) of the input image
            pad       - designates the number of pixels to be padded to the perimeter
                        of the input image.  May be passed as UNDEFINED, in which
                        case the specific padding required by the rotated grids
                        will be computed and returned in ROTGRIDS.
            start_dir_angle - angle from which rotations are to start
            ndirs     - number of rotations to compute (within a semicircle)
            grid_w    - width of the grid in pixels to be rotated
            grid_h    - height of the grid in pixels to be rotated
            relative2 - designates whether pixel offsets whould be computed
                        relative to the ORIGIN or the CENTER of the grid
      **************************************************************************/
         RotGrids(const int iw, const int ipad,
                     const double start_dir_angle, const int ndirs,
                     const int grid_w_, const int grid_h_, const GridOffsetMode relative2_
                  );
      };

      struct Contour {
         int ncontour = 0;
         std::vector<int> xV;
         std::vector<int> yV;
         std::vector<int> exV;
         std::vector<int> eyV;

         Contour()
         {

         }

         explicit Contour(int size)
         {
            Reserve(size);
         }

         void Reserve(int size_) {
            xV.reserve(static_cast<std::size_t>(size_));
            yV.reserve(static_cast<std::size_t>(size_));
            exV.reserve(static_cast<std::size_t>(size_));
            eyV.reserve(static_cast<std::size_t>(size_));
         }

         void AddPoint(int x, int y, int ex, int ey) {
            xV.push_back(x);
            yV.push_back(y);
            exV.push_back(ex);
            eyV.push_back(ey);
            ++ncontour;
         }
      };

      struct FeaturePoint {
         int x;
         int y;
         int ex;
         int ey;
      };

      /*************************************************************************
       DetectMinutiae - Takes a grayscale fingerprint image (of
               arbitrary size), and returns a set of image block maps,
               a binarized image designating ridges from valleys,
               and a list of minutiae (including position, reliability,
               type, direction, neighbors, and ridge counts to neighbors).
               The image maps include a ridge flow directional map,
               a map of low contrast blocks, a map of low ridge flow blocks.
               and a map of high-curvature blocks.

         Input:
            idata     - input 8-bit grayscale fingerprint image data
            iw        - width (in pixels) of the image
            ih        - height (in pixels) of the image
            lqmParams - parameters and thresholds for controlling OpenLQM

         Output:
            ominutiae - resulting list of minutiae
            odmap     - resulting Direction Map
                        {invalid (-1) or valid ridge directions}
            olcmap    - resulting Low Contrast Map
                        {low contrast (TRUE), high contrast (FALSE)}
            olfmap    - resulting Low Ridge Flow Map
                        {low ridge flow (TRUE), high ridge flow (FALSE)}
            ohcmap    - resulting High Curvature Map
                        {high curvature (TRUE), low curvature (FALSE)}
            omw       - width (in blocks) of image maps
            omh       - height (in blocks) of image maps
            obdata    - resulting binarized image
                        {0 = black pixel (ridge) and 255 = white pixel (valley)}
            obw       - width (in pixels) of the binary image
            obh       - height (in pixels) of the binary image
      **************************************************************************/
      void DetectMinutiae(
         Minutiae& ominutiae,
         std::vector<int>& odmap, std::vector<int>& olcmap, std::vector<int>& olfmap, std::vector<int>& ohcmap,
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
      );

      /*************************************************************************
      **************************************************************************
      FindMinutiae - Takes a binary image and its associated
               Direction and Low Flow Maps and scans each image block
               with valid direction for minutia points.  Minutia points
               detected in LOW FLOW blocks are set with lower reliability.

         Input:
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            direction_map  - map of image blocks containing directional ridge flow
            low_flow_map   - map of image blocks flagged as LOW RIDGE FLOW
            high_curve_map - map of image blocks flagged as HIGH CURVATURE
            mw        - width (in blocks) of the maps
            mh        - height (in blocks) of the maps
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae   - reference to a list of detected minutia structures
      **************************************************************************/
      void FindMinutiae(Minutiae& minutiae,
                  unsigned char *bdata, const int iw, const int ih,
                  int *direction_map, int *low_flow_map, int *high_curve_map,
                  const int mw, const int mh,
                  const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      ScanForMinutiaeHorizontally - Scans an entire binary image
                     horizontally, detecting potential minutiae points.
                     Minutia detected via the horizontal scan process are
                     by nature vertically oriented (orthogonal to the scan).

         Input:
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            pdirection_map  - pixelized Direction Map
            plow_flow_map   - pixelized Low Ridge Flow Map
            phigh_curve_map - pixelized High Curvature Map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae   - points to a list of detected minutia structures
      **************************************************************************/
      void ScanForMinutiaeHorizontally(Minutiae& minutiae,
                     unsigned char *bdata, const int iw, const int ih,
                     int *pdirection_map, int *plow_flow_map, int *phigh_curve_map,
                     const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      MatchFirstPair - Determines which of the feature_patterns[] have their
               first pixel pair match the specified pixel pair.

         Input:
            p1 - first pixel value of pair
            p2 - second pixel value of pair
         Output:
            possible - list of matching feature_patterns[] indices
            nposs    - number of matches
         Return Code:
            nposs    - number of matches
      **************************************************************************/
      int MatchFirstPair(unsigned char p1, unsigned char p2, int *possible, int *nposs);

      /*************************************************************************
      **************************************************************************
      MatchSecondPair - Determines which of the passed feature_patterns[] have
               their second pixel pair match the specified pixel pair.

         Input:
            p1 - first pixel value of pair
            p2 - second pixel value of pair
            possible - list of potentially-matching feature_patterns[] indices
            nposs    - number of potential matches
         Output:
            possible - list of matching feature_patterns[] indices
            nposs    - number of matches
         Return Code:
            nposs    - number of matches
      **************************************************************************/
      int MatchSecondPair(unsigned char p1, unsigned char p2, int *possible, int *nposs);

      /*************************************************************************
      **************************************************************************
      SkipRepeatedHorizontalPair - Takes the location of two pixel in
               adjacent pixel rows within an image region and skips
               rightward until the either the pixel pair no longer repeats
               itself or the image region is exhausted.

         Input:
            cx    - current x-coord of starting pixel pair
            ex    - right edge of the image region
            p1ptr - pointer to current top pixel in pair
            p2ptr - pointer to current bottom pixel in pair
         Output:
            cx    - x-coord of where rightward skip terminated
            p1ptr - points to top pixel where rightward skip terminated
            p2ptr - points to bottom pixel where rightward skip terminated
      **************************************************************************/
      void SkipRepeatedHorizontalPair(int *cx, const int ex,
                        unsigned char **p1ptr, unsigned char **p2ptr);

      /*************************************************************************
      **************************************************************************
      SkipRepeatedVerticalPair - Takes the location of two pixel in
               adjacent pixel columns within an image region and skips
               downward until the either the pixel pair no longer repeats
               itself or the image region is exhausted.

         Input:
            cy    - current y-coord of starting pixel pair
            ey    - bottom of the image region
            p1ptr - pointer to current left pixel in pair
            p2ptr - pointer to current right pixel in pair
            iw    - width (in pixels) of image
         Output:
            cy    - y-coord of where downward skip terminated
            p1ptr - points to left pixel where downward skip terminated
            p2ptr - points to right pixel where donward skip terminated
      **************************************************************************/
      void SkipRepeatedVerticalPair(int *cy, const int ey,
                        unsigned char **p1ptr, unsigned char **p2ptr,
                        const int iw);

      /*************************************************************************
      **************************************************************************
      MatchThirdPair - Determines which of the passed feature_patterns[] have
               their third pixel pair match the specified pixel pair.

         Input:
            p1 - first pixel value of pair
            p2 - second pixel value of pair
            possible - list of potentially-matching feature_patterns[] indices
            nposs    - number of potential matches
         Output:
            possible - list of matching feature_patterns[] indices
            nposs    - number of matches
         Return Code:
            nposs    - number of matches
      **************************************************************************/
      int MatchThirdPair(unsigned char p1, unsigned char p2, int *possible, int *nposs);

      /*************************************************************************
      **************************************************************************
      ProcessHorizontalScanMinutia - Takes a minutia point that was
                     detected via the horizontal scan process and
                     adjusts its location (if necessary), determines its
                     direction, and (if it is not already in the minutiae
                     list) adds it to the list.  These minutia are by nature
                     vertical in orientation (orthogonal to the scan).

         Input:
            cx        - x-pixel coord where 3rd pattern pair of mintuia was detected
            cy        - y-pixel coord where 3rd pattern pair of mintuia was detected
            y2        - y-pixel coord where 2nd pattern pair of mintuia was detected
            feature_id - type of minutia (ex. index into feature_patterns[] list)
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            pdirection_map  - pixelized Direction Map
            plow_flow_map   - pixelized Low Ridge Flow Map
            phigh_curve_map - pixelized High Curvature Map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae   - points to a list of detected minutia structures
         Return Code:
            Zero      - successful completion
            LQM_IGNORE - minutia is to be ignored
      **************************************************************************/
      int ProcessHorizontalScanMinutia(Minutiae& minutiae,
                     const int cx, const int cy,
                     const int x2, const int feature_id,
                     unsigned char *bdata, const int iw, const int ih,
                     int *pdirection_map, int *plow_flow_map, int *phigh_curve_map,
                     const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      ProcessVerticalScanMinutia - Takes a minutia point that was
                     detected in via the vertical scan process and
                     adjusts its location (if necessary), determines its
                     direction, and (if it is not already in the minutiae
                     list) adds it to the list.  These minutia are by nature
                     horizontal in orientation (orthogonal to the scan).

         Input:
            cx        - x-pixel coord where 3rd pattern pair of mintuia was detected
            cy        - y-pixel coord where 3rd pattern pair of mintuia was detected
            x2        - x-pixel coord where 2nd pattern pair of mintuia was detected
            feature_id - type of minutia (ex. index into feature_patterns[] list)
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            pdirection_map  - pixelized Direction Map
            plow_flow_map   - pixelized Low Ridge Flow Map
            phigh_curve_map - pixelized High Curvature Map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae   - points to a list of detected minutia structures
         Return Code:
            Zero      - successful completion
            IGNORE    - minutia is to be ignored
      **************************************************************************/
      int ProcessVerticalScanMinutia(Minutiae& minutiae,
                     const int cx, const int cy,
                     const int y2, const int feature_id,
                     unsigned char *bdata, const int iw, const int ih,
                     int *pdirection_map, int *plow_flow_map, int *phigh_curve_map,
                     const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      AdjustHighCurvatureMinutia - Takes an initial minutia point
               in a high-curvature area and adjusts its location and
               direction.  First, it walks and extracts the contour
               of the detected feature looking for and processing any loop
               discovered along the way.  Once the contour is extracted,
               the point of highest-curvature is determined and used to
               adjust the location of the minutia point.  The angle of
               the line perpendicular to the tangent on the high-curvature
               contour at the minutia point is used as the mintutia's
               direction.

         Input:
            x_loc     - starting x-pixel coord of feature (interior to feature)
            y_loc     - starting y-pixel coord of feature (interior to feature)
            x_edge    - x-pixel coord of corresponding edge pixel
                        (exterior to feature)
            y_edge    - y-pixel coord of corresponding edge pixel
                        (exterior to feature)
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            plow_flow_map - pixelized Low Ridge Flow Map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            oidir     - direction of adjusted minutia point
            ox_loc    - adjusted x-pixel coord of feature
            oy_loc    - adjusted y-pixel coord of feature
            ox_edge   - adjusted x-pixel coord of corresponding edge pixel
            oy_edge   - adjusted y-pixel coord of corresponding edge pixel
            minutiae   - points to a list of detected minutia structures
         Return Code:
            Zero      - minutia point processed successfully
            LQM_IGNORE - minutia point is to be ignored
      **************************************************************************/
      int AdjustHighCurvatureMinutia(int *oidir, int *ox_loc, int *oy_loc,
                  int *ox_edge, int *oy_edge,
                  const int x_loc, const int y_loc,
                  const int x_edge, const int y_edge,
                  unsigned char *bdata, const int iw, const int ih,
                  int *plow_flow_map, Minutiae& minutiae, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      GetHighCurvatureContour - Takes the pixel coordinate of a detected
               minutia feature point and its corresponding/adjacent edge
               pixel and attempts to extract a contour of specified length
               of the feature's edge.  The contour is extracted by walking
               the feature's edge a specified number of steps clockwise and
               then counter-clockwise. If a loop is detected while
               extracting the contour, the contour of the loop is returned
               with a return code of (LOOP_FOUND).  If the process fails
               to extract a contour of total specified length, then
               the returned contour length is set to Zero, NO allocated
               memory is returned in this case, and the return code is set
               to Zero.  An alternative implementation would be to return
               the incomplete contour with a return code of (INCOMPLETE).
               For now, NO allocated contour is returned in this case.

         Input:
            half_contour - half the length of the extracted contour
                           (full-length non-loop contour = (half_contourX2)+1)
            x_loc  - starting x-pixel coord of feature (interior to feature)
            y_loc  - starting y-pixel coord of feature (interior to feature)
            x_edge - x-pixel coord of corresponding edge pixel
                     (exterior to feature)
            y_edge - y-pixel coord of corresponding edge pixel
                        (exterior to feature)
            bdata  - binary image data (0==while & 1==black)
            iw     - width (in pixels) of image
            ih     - height (in pixels) of image
         Output:
            ocontour    - Contour object containing point count and lists of interior and exterior coords relative to feature
         Return Code:
            Zero       - resulting contour was successfully extracted or is empty
            LQM_LOOP_FOUND - resulting contour forms a complete loop
      **************************************************************************/
      int GetHighCurvatureContour(Contour& ocontour,
                     const int half_contour,
                     const int x_loc, const int y_loc,
                     const int x_edge, const int y_edge,
                     unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      TraceContour - Takes the pixel coordinate of a detected minutia
               feature point and its corresponding/adjacent edge pixel
               and extracts a contour (up to a specified maximum length)
               of the feature's edge in either a clockwise or counter-
               clockwise direction.  A second point is specified, such that
               if this point is encounted while extracting the contour,
               it is to be assumed that a loop has been found and a code
               of (LOOP_FOUND) is returned with the contour. By independently
               specifying this point, successive calls can be made to
               this routine from the same starting point, and loops across
               successive calls can be detected.

         Input:
            max_len - maximum length of contour to be extracted
            x_loop  - x-pixel coord of point, if encountered, triggers LOOP_FOUND
            y_loop  - y-pixel coord of point, if encountered, triggers LOOP_FOUND
            x_loc   - starting x-pixel coord of feature (interior to feature)
            y_loc   - starting y-pixel coord of feature (interior to feature)
            x_edge  - x-pixel coord of corresponding edge pixel (exterior to feature)
            y_edge  - y-pixel coord of corresponding edge pixel (exterior to feature)
            scan_clock - direction in which neighboring pixels are to be scanned
                     for the next contour pixel
            bdata  - binary image data (0==while & 1==black)
            iw     - width (in pixels) of image
            ih     - height (in pixels) of image
         Output:
            ocontour    - Contour object containing point count and lists of interior and exterior coords relative to feature
         Return Code:
            Zero       - resulting contour was successfully allocated and extracted
            LQM_LOOP_FOUND - resulting contour forms a complete loop
            LQM_IGNORE     - trace is not possible due to state of inputs
      **************************************************************************/
      int TraceContour(Contour& ocontour,
                        const int max_len, const int x_loop, const int y_loop,
                        const int x_loc, const int y_loc,
                        const int x_edge, const int y_edge,
                        const int scan_clock,
                        unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      NextContourPixel - Takes a pixel coordinate of a point determined
               to be on the interior edge of a feature (ridge or valley-
               ending), and attempts to locate a neighboring pixel on the
               feature's contour.  Neighbors of the current feature pixel
               are searched in a specified direction (clockwise or counter-
               clockwise) and the first pair of adjacent/neigboring pixels
               found with the first pixel having the color of the feature
               and the second the opposite color are returned as the next
               point on the contour.  One exception happens when the new
               point is on an "exposed" corner.

         Input:
            cur_x_loc  - x-pixel coord of current point on feature's
                        interior contour
            cur_y_loc  - y-pixel coord of current point on feature's
                        interior contour
            cur_x_edge - x-pixel coord of corresponding edge pixel
                        (exterior to feature)
            cur_y_edge - y-pixel coord of corresponding edge pixel
                        (exterior to feature)
            scan_clock - direction in which neighboring pixels are to be scanned
                        for the next contour pixel
            bdata      - binary image data (0==while & 1==black)
            iw         - width (in pixels) of image
            ih         - height (in pixels) of image
         Output:
            next_x_loc  - x-pixel coord of next point on feature's interior contour
            next_y_loc  - y-pixel coord of next point on feature's interior contour
            next_x_edge - x-pixel coord of corresponding edge (exterior to feature)
            next_y_edge - y-pixel coord of corresponding edge (exterior to feature)
         Return Code:
            LQM_TRUE  - next contour point found and returned
            LQM_FALSE - next contour point NOT found
      **************************************************************************
      **************************************************************************/
      int NextContourPixel(int& next_x_loc, int& next_y_loc,
                     int& next_x_edge, int& next_y_edge,
                     const int cur_x_loc, const int cur_y_loc,
                     const int cur_x_edge, const int cur_y_edge,
                     const int scan_clock,
                     unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      StartScanNeighbor - Takes a two pixel coordinates that are either
               aligned north-to-south or east-to-west, and returns the
               position the second pixel is in realtionship to the first.
               The positions returned are based on 8-connectedness.
               NOTE, this routine does NOT account for diagonal positions.

         Input:
            x_prev - x-coord of first point
            y_prev - y-coord of first point
            x_next - x-coord of second point
            y_next - y-coord of second point
         Return Code:
            LQM_NORTH - second pixel above first
            LQM_SOUTH - second pixel below first
            LQM_EAST  - second pixel right of first
            LQM_WEST  - second pixel left of first
      **************************************************************************/
      int StartScanNeighbor(const int x_prev, const int y_prev, const int x_next, const int y_next);

      /*************************************************************************
      **************************************************************************
      NextScanNeighbor - Advances the given 8-connected neighbor index
               on location in the specifiec direction (clockwise or
               counter-clockwise).

         Input:
            nbr_i      - current 8-connected neighbor index
            scan_clock - direction in which the neighbor index is to be advanced
         Return Code:
            Next neighbor - 8-connected index of next neighbor
      **************************************************************************/
      int NextScanNeighbor(const int nbr_i, const int scan_clock);

      /*************************************************************************
      **************************************************************************
      IsLoopClockwise - Takes a feature's contour points and determines if
               the points are ordered clockwise or counter-clockwise about
               the feature.  The routine also requires a default return
               value be specified in the case the the routine is not able
               to definitively determine the contour's order.  This allows
               the default response to be application-specific.

         Input:
            contour_x   - x-coord list for feature's contour points
            contour_y   - y-coord list for feature's contour points
            ncontour    - number of points in contour
            default_ret - default return code (used when we can't tell the order)
         Return Code:
            LQM_TRUE      - contour determined to be ordered clockwise
            LQM_FALSE     - contour determined to be ordered counter-clockwise
            Default   - could not determine the order of the contour
      **************************************************************************/
      int IsLoopClockwise(const int *contour_x, const int *contour_y, const int ncontour, const int default_ret);

      /*************************************************************************
      **************************************************************************
      ChainCodeLoop - Converts a feature's contour points into an
               8-connected chain code vector.  This encoding represents
               the direction taken between each adjacent point in the
               contour.  Chain codes may be used for many purposes, such
               as computing the perimeter or area of an object, and they
               may be used in object detection and recognition.

         Input:
            contour_x - x-coord list for feature's contour points
            contour_y - y-coord list for feature's contour points
            ncontour  - number of points in contour
         Output:
            ochain    - resulting vector of chain codes
      **************************************************************************/
      void ChainCodeLoop(std::vector<int>& ochain, const int *contour_x, const int *contour_y, const int ncontour);

      /*************************************************************************
      **************************************************************************
      IsChainClockwise - Takes an 8-connected chain code vector and
               determines if the codes are ordered clockwise or
               counter-clockwise.
               The routine also requires a default return value be
               specified in the case the the routine is not able to
               definitively determine the chains direction.  This allows
               the default response to be application-specific.

         Input:
            chain       - chain code vector
            default_ret - default return code (used when we can't tell the order)
         Return Code:
            LQM_TRUE      - chain determined to be ordered clockwise
            LQM_FALSE     - chain determined to be ordered counter-clockwise
            Default       - could not determine the order of the chain
      **************************************************************************/
      int IsChainClockwise(const std::vector<int>& chain, const int default_ret);

      /*************************************************************************
      **************************************************************************
      ProcessMinutiaeLoop - Takes a contour list that has been determined to form
               a complete loop, and processes it. If the loop is sufficiently
               large and elongated, then two minutia points are calculated
               along the loop's longest aspect axis.  If it is determined
               that the loop does not contain minutiae, it is filled in the
               binary image.

         Input:
            contour_x  - x-coord list for loop's contour points
            contour_y  - y-coord list for loop's contour points
            contour_ex - x-coord list for loop's edge points
            contour_ey - y-coord list for loop's edge points
            ncontour   - number of points in contour
            bdata      - binary image data (0==while & 1==black)
            iw         - width (in pixels) of image
            ih         - height (in pixels) of image
            plow_flow_map  - pixelized Low Ridge Flow Map
            lqmParams  - parameters and thresholds for controlling LQM
         Output:
            minutiae    - points to a list of detected minutia structures
            OR
            bdata      - binary image data with loop filled
      **************************************************************************/
      void ProcessMinutiaeLoop(Minutiae& minutiae,
                  const Contour& contour,
                  unsigned char *bdata, const int iw, const int ih,
                  int *plow_flow_map, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      GetLoopAspect - Takes a contour list (determined to form a complete
               loop) and measures the loop's aspect (the largest and smallest
               distances across the loop) and returns the points on the
               loop where these distances occur.

         Input:
            contour_x - x-coord list for loop's contour points
            contour_y - y-coord list for loop's contour points
            ncontour  - number of points in contour
         Output:
            omin_fr   - contour point index where minimum aspect occurs
            omin_to   - opposite contour point index where minimum aspect occurs
            omin_dist - the minimum distance across the loop
            omax_fr   - contour point index where maximum aspect occurs
            omax_to   - contour point index where maximum aspect occurs
            omax_dist - the maximum distance across the loop
      **************************************************************************/
      void GetLoopAspect(int *omin_fr, int *omin_to, double *omin_dist,
                  int *omax_fr, int *omax_to, double *omax_dist,
                  const Contour& contour);

      /*************************************************************************
      **************************************************************************
      SquaredDistance - Takes two coordinate points and computes the
                        squared distance between the two points.

         Input:
            x1  - x-coord of first point
            y1  - y-coord of first point
            x2  - x-coord of second point
            y2  - y-coord of second point
         Return Code:
            Distance - computed squared distance
      **************************************************************************/
      template <typename T>
      double SquaredDistance(const T x1, const T y1, const T x2, const T y2)
      {
         double dx, dy, dist;

         /* Compute delta x between points. */
         dx = static_cast<double>(x1 - x2);
         /* Compute delta y between points. */
         dy = static_cast<double>(y1 - y2);
         /* Compute the squared distance between points. */
         dist = (dx*dx) + (dy*dy);

         /* Return the squared distance. */
         return dist;
      }

      /*************************************************************************
      **************************************************************************
      LineToDirection - Takes two coordinate points and computes the
               directon (on a full circle) in which the first points
               to the second.

         Input:
            fx         - x-coord of first point (pointing from)
            fy         - y-coord of first point (pointing from)
            tx         - x-coord of second point (pointing to)
            ty         - y-coord of second point (pointing to)
            ndirs      - number of IMAP directions (in semicircle)
         Return Code:
            Direction  - determined direction on a "full" circle
      **************************************************************************/
      int LineToDirection(const int fx, const int fy, const int tx, const int ty, const int ndirs);

      /*************************************************************************
      **************************************************************************
      AngleToLine - Takes two coordinate points and computes the angle
               to the line formed by the two points.

         Input:
            fx         - x-coord of first point
            fy         - y-coord of first point
            tx         - x-coord of second point
            ty         - y-coord of second point
         Return Code:
            Angle - angle to the specified line
      **************************************************************************/
      double AngleToLine(const int fx, const int fy, const int tx, const int ty);

      /*************************************************************************
      **************************************************************************
      IsMinutiaAppearing - Given the pixel location of a minutia feature
               and its corresponding adjacent edge pixel, returns whether
               the minutia is appearing or disappearing.  Remeber, that
               "feature" refers to either a ridge or valley-ending.

         Input:
            x_loc      - x-pixel coord of feature (interior to feature)
            y_loc      - y-pixel coord of feature (interior to feature)
            x_edge     - x-pixel coord of corresponding edge pixel
                        (exterior to feature)
            y_edge     - y-pixel coord of corresponding edge pixel
                        (exterior to feature)
         Return Code:
            LQM_APPEARING    - minutia is appearing (TRUE==1)
            LQM_DISAPPEARING - minutia is disappearing (FALSE==0)
      **************************************************************************/
      int IsMinutiaAppearing(const int x_loc, const int y_loc, const int x_edge, const int y_edge);

      /*************************************************************************
      **************************************************************************
      UpdateMinutiae - Takes a detected minutia point and (if it is not
                     determined to already be in the minutiae list) adds it to
                     the list.

         Input:
            minutia   - minutia structure for detected point
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae   - points to a list of detected minutia structures
         Return Code:
            Zero      - minutia added to successfully added to minutiae list
            LQM_IGNORE    - minutia is to be ignored (already in the minutiae list)
      **************************************************************************/
      int UpdateMinutiae(Minutiae& minutiae, Minutia& minutia,
                        unsigned char *bdata, const int iw, const int ih,
                        const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      SearchContour - Walk the contour of a minutia feature starting at a
               specified point on the feature and walking N steps in the
               specified direction (clockwise or counter-clockwise), looking
               for a second specified point.  In this code, "feature" is
               consistently referring to either the black interior edge of
               a ridge-ending or the white interior edge of a valley-ending
               (bifurcation).  The term "edge of the feature" refers to
               neighboring pixels on the "exterior" edge of the feature.
               So "edge" pixels are opposite in color from the interior
               feature pixels.

         Input:
            x_search   - x-pixel coord of point being searched for
            y_search   - y-pixel coord of point being searched for
            search_len - number of step to walk contour in search
            x_loc      - starting x-pixel coord of feature (interior to feature)
            y_loc      - starting y-pixel coord of feature (interior to feature)
            x_edge     - x-pixel coord of corresponding edge pixel
                        (exterior to feature)
            y_edge     - y-pixel coord of corresponding edge pixel
                        (exterior to feature)
            scan_clock - direction in which neighbor pixels are to be scanned
                        (clockwise or counter-clockwise)
            bdata      - binary image data (0==while & 1==black)
            iw         - width (in pixels) of image
            ih         - height (in pixels) of image
         Return Code:
            LQM_NOT_FOUND  - desired pixel not found along N steps of feature's contour
            LQM_FOUND      - desired pixel WAS found along N steps of feature's contour
      **************************************************************************/
      int SearchContour(const int x_search, const int y_search,
                        const int search_len,
                        const int x_loc, const int y_loc,
                        const int x_edge, const int y_edge,
                        const int scan_clock,
                        unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      FillLoop - Takes a contour list that has been determined to form
               a complete loop, and fills the loop accounting for
               complex/concaved shapes.
               NOTE, I tried using a flood-fill in place of this routine,
               but the contour (although 8-connected) is NOT guaranteed to
               be "complete" surrounded (in an 8-connected sense) by pixels
               of opposite color.  Therefore, the flood would occasionally
               escape the loop and corrupt the binary image!

         Input:
            contour_x  - x-coord list for loop's contour points
            contour_y  - y-coord list for loop's contour points
            ncontour   - number of points in contour
            bdata      - binary image data (0==while & 1==black)
            iw         - width (in pixels) of image
         Output:
            bdata      - binary image data with loop filled
      **************************************************************************/
      void FillLoop(const int *contour_x, const int *contour_y,
                     const int ncontour, unsigned char *bdata,
                     const int iw);

      /*************************************************************************
      **************************************************************************
      FillPartialRow - Fills a specified range of contiguous pixels on
               a specified row of an 8-bit pixel image with a specified
               pixel value.  NOTE, the pixel coordinates are assumed to
               be within the image boundaries.

         Input:
            fill_pix - pixel value to fill with (should be on range [0..255]
            frx      - x-pixel coord where fill should begin
            tox      - x-pixel coord where fill should end (inclusive)
            y        - y-pixel coord of current row being filled
            bdata    - 8-bit image data
            iw       - width (in pixels) of image
            ih       - height (in pixels) of image
         Output:
            bdata    - 8-bit image data with partial row filled.
      **************************************************************************/
      void FillPartialRow(const int fill_pix, const int frx, const int tox,
               const int y, unsigned char *bdata, const int iw);

      /*************************************************************************
      **************************************************************************
      MinContourTheta - Takes a contour list and analyzes it locating the
               point at which the contour has highest curvature
               (or minimum interior angle).  The angle of curvature is
               computed by searching a majority of points on the contour.
               At each of these points, a left and right segment (or edge)
               are extended out N number of pixels from the center point
               on the contour.  The angle is formed between the straight line
               connecting the center point to the end point on the left edge
               and the line connecting the center point to the end of the
               right edge.  The point of highest curvature is determined
               by locating the where the minimum of these angles occurs.

         Input:
            angle_edge - length of the left and right edges extending from a
                        common/centered pixel on the contour
            contour_x  - x-coord list for contour points
            contour_y  - y-coord list for contour points
            ncontour   - number of points in contour
         Output:
            omin_i     - index of contour point where minimum occurred
            omin_theta - minimum angle found along the contour
         Return Code:
            Zero       - minimum angle successfully located
            LQM_IGNORE - ignore the contour
      **************************************************************************/
      int MinContourTheta(int *omin_i, double *omin_theta,
                           const int angle_edge,  const int *contour_x,
                           const int *contour_y, const int ncontour);

      /*************************************************************************
      **************************************************************************
      GetLowCurvatureDirection - Converts a bi-direcitonal IMAP direction
               (based on a semi-circle) to a uni-directional value covering
               a full circle based on the scan orientation used to detect
               a minutia feature (horizontal or vertical) and whether the
               detected minutia is appearing or disappearing.

         Input:
            scan_dir   - designates the feature scan orientation
            appearing  - designates the minutia as appearing or disappearing
            imapval    - IMAP block direction
            ndirs      - number of IMAP directions (in semicircle)
         Return Code:
            New direction - bi-directonal integer direction on full circle
      **************************************************************************/
      int GetLowCurvatureDirection(const int scan_dir, const int appearing, const int imapval, const int ndirs);

      /*************************************************************************
      **************************************************************************
      ScanAndUpdateMinutiae - Takes a detected minutia point and (if it is not
                     determined to already be in the minutiae list or the
                     new point is determined to be "more compatible") adds
                     it to the list.

         Input:
            minutia   - minutia structure for detected point
            scan_dir  - orientation of scan when minutia was detected
            dmapval   - directional ridge flow of block minutia is in
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae   - reference to a list of detected minutia structures
         Return Code:
            Zero      - minutia added to successfully added to minutiae list
            LQM_IGNORE - minutia is to be ignored (already in the minutiae list)
      **************************************************************************/
      int ScanAndUpdateMinutiae(Minutiae& minutiae, Minutia& minutia,
                        const int scan_dir, const int dmapval,
                        unsigned char *bdata, const int iw, const int ih,
                        const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      ChooseScanDirection - Determines the orientation (horizontal or
               vertical) in which a block is to be scanned for minutiae.
               The orientation is based on the blocks corresponding IMAP
               direction.

         Input:
            imapval   - Block's IMAP direction
            ndirs     - number of possible IMAP directions (within semicircle)
         Return Code:
            LQM_SCAN_HORIZONTAL - horizontal orientation
            LQM_SCAN_VERTICAL   - vertical orientation
      **************************************************************************/
      int ChooseScanDirection(const int imapval, const int ndirs);

      /*************************************************************************
      **************************************************************************
      ScanForMinutiaeVertically - Scans an entire binary image
                     vertically, detecting potential minutiae points.
                     Minutia detected via the vetical scan process are
                     by nature horizontally oriented (orthogonal to  the scan).

         Input:
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            pdirection_map  - pixelized Direction Map
            plow_flow_map   - pixelized Low Ridge Flow Map
            phigh_curve_map - pixelized High Curvature Map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae   - reference to a list of detected minutia structures
      **************************************************************************/
      void ScanForMinutiaeVertically(Minutiae& minutiae,
                     unsigned char *bdata, const int iw, const int ih,
                     int *pdirection_map, int *plow_flow_map, int *phigh_curve_map,
                     const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      RemoveFalseMinutia - Takes a list of true and false minutiae and
                     attempts to detect and remove the false minutiae based
                     on a series of tests.

         Input:
            minutiae  - list of true and false minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            direction_map  - map of image blocks containing directional ridge flow
            low_flow_map   - map of image blocks flagged as LOW RIDGE FLOW
            high_curve_map - map of image blocks flagged as HIGH CURVATURE
            mw        - width in blocks of the maps
            mh        - height in blocks of the maps
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemoveFalseMinutia(Minutiae& minutiae,
               unsigned char *bdata, const int iw, const int ih,
               int *direction_map, int *low_flow_map, int *high_curve_map,
               const int mw, const int mh, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      SortMinutiaeYX - Takes a list of minutia points and sorts them
                     top-to-bottom and then left-to-right.

         Input:
            minutiae  - list of minutiae
            iw        - width (in pixels) of image
         Output:
            minutiae  - list of sorted minutiae
      **************************************************************************/
      void SortMinutiaeYX(Minutiae& minutiae, const int iw);

      /*************************************************************************
      **************************************************************************
      SortMinutiaeXY - Takes a list of minutia points and sorts them
                     left-to-right and then top-to-bottom.

         Input:
            minutiae  - list of minutiae
            iw        - width (in pixels) of image
         Output:
            minutiae  - list of sorted minutiae
      **************************************************************************/
      void SortMinutiaeXY(Minutiae& minutiae, const int iw);

      /*************************************************************************
      **************************************************************************
      RemoveIslandsAndLakes - Takes a list of true and false minutiae and
                     attempts to detect and remove those false minutiae that
                     are either on a common island (filled with black pixels)
                     or a lake (filled with white pixels).
                     Note that this routine edits the binary image by filling
                     detected lakes or islands.

         Input:
            minutiae  - list of true and false minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemoveIslandsAndLakes(Minutiae& minutiae,
                           unsigned char *bdata, const int iw, const int ih,
                           const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      OnIslandLake - Determines if two minutia points lie on the same loop
                     (island or lake).  If a loop is detected, the contour
                     points of the loop are returned.

         Input:
            minutia1      - first minutia point
            minutia2      - second minutia point
            max_half_loop - maximum size of half the loop circumference searched for
            bdata         - binary image data (0==while & 1==black)
            iw            - width (in pixels) of image
            ih            - height (in pixels) of image
         Output:
            ocontour    - Loop contour object containing point count and lists of interior and exterior coords
         Return Code:
            LQM_IGNORE     - contour could not be traced
            LQM_LOOP_FOUND - minutiae determined to lie on same qualifying loop
            LQM_FALSE      - minutiae determined not to lie on same qualifying loop
      **************************************************************************/
      int OnIslandLake(Contour& ocontour,
                        const Minutia& minutia1, const Minutia& minutia2,
                        const int max_half_loop,
                        unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      RemoveHoles - Removes minutia points on small loops around valleys.

         Input:
            minutiae  - list of true and false minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemoveHoles(Minutiae& minutiae,
                     unsigned char *bdata, const int iw, const int ih,
                     const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      MinutiaOnLoop - Determines if a minutia point lies on a loop (island or lake)
               of specified maximum circumference.

         Input:
            minutiae      - list of true and false minutiae
            max_loop_len  - maximum size of loop searched for
            bdata         - binary image data (0==while & 1==black)
            iw            - width (in pixels) of image
            ih            - height (in pixels) of image
         Return Code:
            LQM_IGNORE     - minutia contour could not be traced
            LQM_LOOP_FOUND - minutia determined to lie on qualifying loop
            LQM_FALSE      - minutia determined not to lie on qualifying loop
      **************************************************************************/
      int MinutiaOnLoop(const Minutia& minutia, const int max_loop_len,
                  unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      RemovePointingInvalidBlock - Removes minutia points that are relatively
                     close in the direction opposite the minutia to a
                     block with INVALID ridge flow.

         Input:
            minutiae  - list of true and false minutiae
            direction_map - map of image blocks containing directional ridge flow
            mw        - width in blocks of the map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemovePointingInvalidBlock(Minutiae& minutiae,
                                 int *direction_map, const int mw,
                                 const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      RemoveNearInvalidBlocks - Removes minutia points from the given list
                     that are sufficiently close to a block with invalid
                     ridge flow or to the edge of the image.

         Input:
            minutiae  - list of true and false minutiae
            direction_map - map of image blocks containing direction ridge flow
            mw        - width in blocks of the map
            mh        - height in blocks of the map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemoveNearInvalidBlocks(Minutiae& minutiae, int *direction_map, const int mw, const int mh, const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      RemoveOrAdjustSideMinutiae - Removes loops or minutia points that
                  are not on complete contours of specified length. If the
                  contour is complete, then the minutia is adjusted based
                  on a minmax analysis of the rotated y-coords of the contour.

         Input:
            minutiae  - list of true and false minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            direction_map - map of image blocks containing directional ridge flow
            mw        - width (in blocks) of the map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemoveOrAdjustSideMinutiae(Minutiae& minutiae,
                     unsigned char *bdata, const int iw, const int ih,
                     int *direction_map, const int mw,
                     const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      GetCenteredContour - Takes the pixel coordinate of a detected
               minutia feature point and its corresponding/adjacent edge
               pixel and attempts to extract a contour of specified length
               of the feature's edge.  The contour is extracted by walking
               the feature's edge a specified number of steps clockwise and
               then counter-clockwise. If a loop is detected while
               extracting the contour, no contour is returned with a return
               code of (LOOP_FOUND).  If the process fails to extract a
               a complete contour, a code of INCOMPLETE is returned.

         Input:
            half_contour - half the length of the extracted contour
                           (full-length non-loop contour = (half_contourX2)+1)
            x_loc  - starting x-pixel coord of feature (interior to feature)
            y_loc  - starting y-pixel coord of feature (interior to feature)
            x_edge - x-pixel coord of corresponding edge pixel
                     (exterior to feature)
            y_edge - y-pixel coord of corresponding edge pixel
                        (exterior to feature)
            bdata  - binary image data (0==while & 1==black)
            iw     - width (in pixels) of image
            ih     - height (in pixels) of image
         Output:
            ocontour    - Contour object containing point count and lists of interior and exterior coords relative to feature
         Return Code:
            Zero       - resulting contour was successfully extracted or is empty
            LQM_LOOP_FOUND - resulting contour forms a complete loop
            LQM_IGNORE     - contour could not be traced due to problem starting
                        conditions
            LQM_INCOMPLETE - resulting contour was not long enough
      **************************************************************************/
      int GetCenteredContour(Contour& ocontour,
                     const int half_contour,
                     const int x_loc, const int y_loc,
                     const int x_edge, const int y_edge,
                     unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      MINMAXS - Takes a list of integers and identifies points of relative
               minima and maxima.  The midpoint of flat plateaus and valleys
               are selected when they are detected.

         Input:
            items     - list of integers to be analyzed
            num       - number of items in the list
         Output:
            ominmax_val   - value of the item at each minima or maxima
            ominmax_type  - identifies a minima as '-1' and maxima as '1'
            ominmax_i     - index of item's position in list
            ominmax_alloc - number of allocated minima and/or maxima
            ominmax_num   - number of detected minima and/or maxima
      **************************************************************************/
      void GetExtrema(std::vector<Extremum>& extrema, const int *items, const int num);

      /*************************************************************************
      **************************************************************************
      RemoveHooks - Takes a list of true and false minutiae and
                     attempts to detect and remove those false minutiae that
                     are on a hook (white or black).

         Input:
            minutiae  - list of true and false minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemoveHooks(Minutiae& minutiae,
                     unsigned char *bdata, const int iw, const int ih,
                     const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      OnHook - Determines if two minutia points lie on a hook on the side
               of a ridge or valley.

         Input:
            minutia1      - first minutia point
            minutia2      - second minutia point
            max_hook_len  - maximum length of contour searched along for a hook
            bdata         - binary image data (0==while & 1==black)
            iw            - width (in pixels) of image
            ih            - height (in pixels) of image
         Return Code:
            LQM_IGNORE     - contour could not be traced
            LQM_HOOK_FOUND - minutiae determined to lie on same qualifying hook
            LQM_FALSE      - minutiae determined not to lie on same qualifying hook
      **************************************************************************/
      int OnHook(const Minutia& minutia1, const Minutia& minutia2,
                  const int max_hook_len,
                  unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      RemoveOverlaps - Takes a list of true and false minutiae and
                     attempts to detect and remove those false minutiae that
                     are on opposite sides of an overlap.  Note that this
                     routine does NOT edit the binary image when overlaps
                     are removed.

         Input:
            minutiae  - list of true and false minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemoveOverlaps(Minutiae& minutiae,
                        unsigned char *bdata, const int iw,
                        const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      RemoveMalformations - Attempts to detect and remove minutia points
               that are "irregularly" shaped.  Irregularity is measured
               by measuring across the interior of the feature at
               two progressive points down the feature's contour.  The
               test is triggered if a pixel of opposite color from the
               feture's type is found.  The ratio of the distances across
               the feature at the two points is computed and if the ratio
               is too large then the minutia is determined to be malformed.
               A cursory test is conducted prior to the general tests in
               the event that the minutia lies in a block with LOW RIDGE
               FLOW.  In this case, the distance across the feature at
               the second progressive contour point is measured and if
               too large, the point is determined to be malformed.

         Input:
            minutiae  - list of true and false minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            low_flow_map   - map of image blocks flagged as LOW RIDGE FLOW
            mw        - width in blocks of the map
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemoveMalformations(Minutiae& minutiae,
                              unsigned char *bdata, const int iw, const int ih,
                              int *low_flow_map, const int mw,
                              const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      RemovePores - Attempts to detect and remove minutia points located on
                        pore-shaped valleys and/or ridges.  Detection for
                        these features are only performed in blocks with
                        LOW RIDGE FLOW or HIGH CURVATURE.

         Input:
            minutiae  - list of true and false minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            direction_map  - map of image blocks containing directional ridge flow
            low_flow_map   - map of image blocks flagged as LOW RIDGE FLOW
            high_curve_map - map of image blocks flagged as HIGH CURVATURE
            mw        - width in blocks of the maps
            lqmParams - parameters and thresholds for controlling LQM
         Output:
            minutiae  - list of pruned minutiae
      **************************************************************************/
      void RemovePores(Minutiae& minutiae,
                        unsigned char *bdata, const int iw, const int ih,
                        int *direction_map, int *low_flow_map,
                        int *high_curve_map, const int mw,
                        const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      SearchInDirection - Takes a specified maximum number of steps in a
                     specified direction looking for the first occurence of
                     a pixel with specified value.  (Once found, adjustments
                     are potentially made to make sure the resulting pixel
                     and its associated edge pixel are 4-connected.)

         Input:
            pix       - value of pixel to be searched for
            strt_x    - x-pixel coord to start search
            strt_y    - y-pixel coord to start search
            delta_x   - increment in x for each step
            delta_y   - increment in y for each step
            maxsteps  - maximum number of steps to conduct search
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
         Output:
            op        - FeaturePoint object containing coordinates of the located pixel and its associated edge pixel
            ox        - x coord of located pixel
            oy        - y coord of located pixel
            oex       - x coord of associated edge pixel
            oey       - y coord of associated edge pixel
         Return Code:
            true      - pixel of specified value found
            false     - pixel of specified value NOT found
      **************************************************************************/
      bool SearchInDirection(FeaturePoint& op, const int pix,
                     const int strt_x, const int strt_y,
                     const double delta_x, const double delta_y, const int maxsteps,
                     unsigned char *bdata, const int iw, const int ih);

      /*************************************************************************
      **************************************************************************
      FixEdgePixelPair - Takes a pair of pixel points with the first
            pixel on a feature and the second adjacent and off the feature,
            determines if the pair neighbor diagonally.  If they do, their
            locations are adjusted so that the resulting pair retains the
            same pixel values, but are neighboring either to the N,S,E or W.
            This routine is needed in order to prepare the pixel pair for
            contour tracing.

         Input:
            p       - reference to feature point object
            bdata   - binary image data (0==while & 1==black)
            iw      - width (in pixels) of image
         Output:
            p       - reference to resulting feature point object
      **************************************************************************/
      void FixEdgePixelPair(FeaturePoint& p, unsigned char *bdata, const int iw);

      /*************************************************************************
      **************************************************************************
      CountMinutiaeRidges - Takes a list of minutiae, and for each one,
                     determines its closest neighbors and counts the number
                     of interveining ridges between the minutia point and
                     each of its neighbors.

         Input:
            minutiae  - list of minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            lqmParams - parameters and thresholds specified by LQM
         Output:
            minutiae  - list of minutiae augmented with neighbors and ridge counts
      **************************************************************************/
      void CountMinutiaeRidges(Minutiae& minutiae,
         unsigned char *bdata, const int iw, const int ih,
         const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      CountMinutiaRidges - Takes a minutia, and determines its closest
                     neighbors and counts the number of interveining ridges
                     between the minutia point and each of its neighbors.

         Input:
            minutia   - input minutia
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            lqmParams - parameters and thresholds specified by LQM
         Output:
            minutiae  - minutia augmented with neighbors and ridge counts
      **************************************************************************/
      void CountMinutiaRidges(const int first, Minutiae& minutiae,
                              unsigned char *bdata, const int iw, const int ih,
                              const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      FindNeighbors - Takes a primary minutia and a list of all minutiae
                  and locates a specified maximum number of closest neighbors
                  to the primary point.  Neighbors are searched, starting
                  in the same pixel column, below, the primary point and then
                  along consecutive and complete pixel columns in the image
                  to the right of the primary point.

         Input:
            max_nbrs - maximum number of closest neighbors to be returned
            first    - index of the primary minutia point
            minutiae - list of minutiae
         Output:
            onbr_list - list of detected closest neighbors
      **************************************************************************/
      void FindNeighbors(std::vector<int>& onbr_list, const int max_nbrs, const int first, Minutiae& minutiae);

      /*************************************************************************
      **************************************************************************
      UpdateNeighborDists - Takes the current list of neighbors along with a
                  primary minutia and a potential new neighbor, and
                  determines if the new neighbor is sufficiently close
                  to be added to the list of nearest neighbors.  If added,
                  it is placed in the list in its proper order based on
                  squared distance to the primary point.

         Input:
            nbr_list - current list of nearest neighbor minutia indices
            nbr_sqr_dists - corresponding squared euclidean distance of each
                     neighbor to the primary minutia point
            max_nbrs - maximum number of closest neighbors to be returned
            first    - index of the primary minutia point
            second   - index of the secondary (new neighbor) point
            minutiae - list of minutiae
         Output:
            nbr_list - updated list of nearest neighbor indices
            nbr_sqr_dists - updated list of nearest neighbor distances
      **************************************************************************/
      void UpdateNeighborDists(std::vector<int>& nbr_list, std::vector<double>& nbr_sqr_dists,
         const int max_nbrs, const int first, const int second, Minutiae& minutiae);

      /*************************************************************************
      **************************************************************************
      FindIncreasingPositionDouble - Takes a double value and a list of doubles and
                  determines where in the list the double may be inserted,
                  preserving the increasing sorted order of the list.

         Input:
            val  - value to be inserted into the list
            list - list of double in increasing sorted order
         Return Code:
            Zero or Positive - insertion position in the list
      **************************************************************************/
      int FindIncreasingPositionDouble(const double val, std::vector<double>& list);

      /*************************************************************************
      **************************************************************************
      InsertNeighbor - Takes a minutia index and its squared distance to a
                  primary minutia point, and inserts them in the specified
                  position of their respective lists, shifting previously
                  stored values down and off the lists as necessary.

         Input:
            pos       - postions where values are to be inserted in lists
            nbr_index - index of minutia being inserted
            nbr_dist2 - squared distance of minutia to its primary point
            nbr_list  - current list of nearest neighbor minutia indices
            nbr_sqr_dists - corresponding squared euclidean distance of each
                        neighbor to the primary minutia point
            max_nbrs  - maximum number of closest neighbors to be returned
         Output:
            nbr_list - updated list of nearest neighbor indices
            nbr_sqr_dists - updated list of nearest neighbor distances
      **************************************************************************/
      void InsertNeighbor(const int pos, const int nbr_index, const double nbr_dist2,
         std::vector<int>& nbr_list, std::vector<double>& nbr_sqr_dists, const int max_nbrs);

      /*************************************************************************
      **************************************************************************
      SortNeighbors - Takes a list of primary minutia and its neighboring
                  minutia indices and sorts the neighbors based on their
                  position relative to the primary minutia point.  Neighbors
                  are sorted starting vertical to the primary point and
                  proceeding clockwise.

         Input:
            nbr_list - list of neighboring minutia indices
            first    - the index of the primary minutia point
            minutiae - list of minutiae
         Output:
            nbr_list - neighboring minutia indices in sorted order
      **************************************************************************/
      void SortNeighbors(std::vector<int>& nbr_list, const int first, Minutiae& minutiae);

      /*************************************************************************
      **************************************************************************
      RidgeCount - Takes a pair of minutiae, and counts the number of
                  ridges crossed along the linear trajectory connecting
                  the 2 points in the image.

         Input:
            first     - index of primary minutia
            second    - index of secondary (neighbor) minutia
            minutiae  - list of minutiae
            bdata     - binary image data (0==while & 1==black)
            iw        - width (in pixels) of image
            ih        - height (in pixels) of image
            lqmParams - parameters and thresholds specified by LQM
         Return Code:
            Zero or Positive - number of ridges counted
      **************************************************************************/
      int RidgeCount(const int first, const int second, Minutiae& minutiae,
         unsigned char *bdata, const int iw, const int ih,
         const LQMParams& lqmParams);

      /*************************************************************************
      **************************************************************************
      FindTransition - Takes a pixel trajectory and a starting index, and
                  searches forward along the trajectory until the specified
                  adjacent pixel pair is found, returning the index where
                  the pair was found (the index of the second pixel).

         Input:
            iptr  - pointer to starting pixel index into trajectory
            pix1  - first pixel value in transition pair
            pix2  - second pixel value in transition pair
            points - list of pixel coords of line trajectory
            bdata - binary image data (0==while & 1==black)
            iw    - width (in pixels) of image
         Output:
            iptr  - points to location where 2nd pixel in pair is found
         Return Code:
            true  - pixel pair transition found
            false - pixel pair transition not found
      **************************************************************************/
      bool FindTransition(int *iptr, const int pix1, const int pix2,
         const std::vector<std::pair<int, int>>& points,
         unsigned char *bdata, const int iw);

      /*************************************************************************
      **************************************************************************
      ValidateRidgeCrossing - Takes a pair of points, a ridge start
                  transition and a ridge end transition, and walks the
                  ridge contour from thre ridge end points a specified
                  number of steps, looking for the ridge start point.
                  If found, then transitions determined not to be a valid
                  ridge crossing.

         Input:
            ridge_start - index into line trajectory of ridge start transition
            ridge_end   - index into line trajectory of ridge end transition
            points      - list of pixel coords of line trajectory
            bdata       - binary image data (0==while & 1==black)
            iw          - width (in pixels) of image
            ih          - height (in pixels) of image
            max_ridge_steps  - number of steps taken in search in both
                              scan directions
         Return Code:
            true        - ridge crossing VALID
            false       - ridge corssing INVALID
      **************************************************************************/
      bool ValidateRidgeCrossing(const int ridge_start, const int ridge_end,
         const std::vector<std::pair<int, int>>& points,
         unsigned char *bdata, const int iw, const int ih,
         const int max_ridge_steps);

      /*************************************************************************
      **************************************************************************
      RemoveDuplicateMinutiae - Takes a list of minutiae sorted in some adjacent order
                  and detects and removes redundant minutia that have the
                  same exact pixel coordinate locations (even if other
                  attributes may differ).

         Input:
            minutiae - list of sorted minutiae
         Output:
            minutiae - list of sorted minutiae with duplicates removed
      **************************************************************************/
      void RemoveDuplicateMinutiae(Minutiae& minutiae);

      int FilterAndOutputMinutiae(MinutiaOut* pOut, int maxMinutiae, Minutiae& minutiae, int inputResolution);

      bool FindCloseMinutiae(const Minutiae& minutiae, const Minutia& testMinutia, int pixelDistance);

      long AngleToDegrees(int angle);

      // Given a point, an angle (ccw from right), and a distance, returns a point that distance and angle from the point
      std::pair<int, int> GetDirectedPoint(int x, int y, double angle, double length);

      void FilterMinutiaeByQuality(const unsigned char* pLocQ, int mapWidth, int mapHeight, int inputWidth, int inputHeight, MinutiaOut* pMinutiae, int minutiaCount, std::vector<MinutiaOut>& sortedMinutiae);

      void CalculateMinutiaCounts(const std::vector<MinutiaOut>& minutiae, std::array<int, 3>& minutiaCounts);
   }
}

namespace OpenLQM {
    /* Constants (C) for defining 4 DFT frequencies, where  */
    /* frequency is defined as C*(PI_FACTOR).  PI_FACTOR    */
    /* regulates the period of the function in x, so:       */
    /*      1 = one period in range X.                      */
    /*      2 = twice the frequency in range X.             */
    /*      3 = three times the frequency in reange X.      */
    /*      4 = four times the frequency in ranage X.       */
    inline double DFT_COEFFICIENTS[NUM_DFT_WAVES] = { 1,2,3,4 };

    /* Global array of feature pixel pairs. */
    inline OpenLQM::Core::FeaturePattern FEATURE_PATTERNS[]=
                        {{LQM_RIDGE_ENDING,  /* a. Ridge Ending (appearing) */
                            LQM_APPEARING,
                            {0,0},
                            {0,1},
                            {0,0}},

                            {LQM_RIDGE_ENDING,  /* b. Ridge Ending (disappearing) */
                            LQM_DISAPPEARING,
                            {0,0},
                            {1,0},
                            {0,0}},

                            {LQM_BIFURCATION,   /* c. Bifurcation (disappearing) */
                            LQM_DISAPPEARING,
                            {1,1},
                            {0,1},
                            {1,1}},

                            {LQM_BIFURCATION,   /* d. Bifurcation (appearing) */
                            LQM_APPEARING,
                            {1,1},
                            {1,0},
                            {1,1}},

                            {LQM_BIFURCATION,   /* e. Bifurcation (disappearing) */
                            LQM_DISAPPEARING,
                            {1,0},
                            {0,1},
                            {1,1}},

                            {LQM_BIFURCATION,   /* f. Bifurcation (disappearing) */
                            LQM_DISAPPEARING,
                            {1,1},
                            {0,1},
                            {1,0}},

                            {LQM_BIFURCATION,   /* g. Bifurcation (appearing) */
                            LQM_APPEARING,
                            {1,1},
                            {1,0},
                            {0,1}},

                            {LQM_BIFURCATION,   /* h. Bifurcation (appearing) */
                            LQM_APPEARING,
                            {0,1},
                            {1,0},
                            {1,1}},

                            {LQM_BIFURCATION,   /* i. Bifurcation (disappearing) */
                            LQM_DISAPPEARING,
                            {1,0},
                            {0,1},
                            {1,0}},

                            {LQM_BIFURCATION,   /* j. Bifurcation (appearing) */
                            LQM_APPEARING,
                            {0,1},
                            {1,0},
                            {0,1}}};

    /* Variables for conducting 8-connected neighbor analyses. */
    /* Pixel neighbor offsets:  0  1  2  3  4  5  6  7  */     /* 7 0 1 */
    constexpr inline int NEIGH8_DX[]{  0, 1, 1, 1, 0,-1,-1,-1 };      /* 6 C 2 */
    constexpr inline int NEIGH8_DY[]{ -1,-1, 0, 1, 1, 1, 0,-1 };      /* 5 4 3 */

    /* The chain code lookup matrix for 8-connected neighbors. */
    constexpr inline int NEIGH8_DIM = 3;
    constexpr inline int NEIGH8_CHAINCODES[] = {
        3, 2, 1,
        4,-1, 0,
        5, 6, 7
    };

    namespace Presets {
        namespace MinutiaeReliability {
            constexpr inline double MEDIUM = 0.50;
            constexpr inline double HIGH = 0.99;
            constexpr inline double DEFAULT = HIGH;
        }

        constexpr inline OpenLQM::Core::LQMParams LQM_PARAMS = {
            128,
            1,
            4,
            24,
            10,
            16,
            LQM_PI/2.0,
            3,
            0.2,
            3,
            7,
            7,
            5,
            5,
            2,
            10,
            5,
            4,
            100000.0,
            3.8,
            50000000.0,
            2,
            0.7,
            0.75,
            7,
            9,
            0,
            3,
            10,
            LQM_PI/3.0,
            14,
            20,
            1.0,
            2.25,
            0,
            0,
            0,
            2,
            0.0,
            0.0,
            0.0,
            0.0,
            16,
            30,
            30,
            2,
            15,
            7,
            2,
            7,
            8,
            6,
            10,
            20,
            2.0,
            20,
            3,
            12,
            10,
            8,
            0.5,
            2.25,
            5,
            10
        };
    }

    namespace Constants {
        /***** MINUTIAE DETECTION CONSTANTS *****/
        namespace Minutia {
            /* There are 10 unique feature patterns with ID = [0..9] , */
            /* so set LOOP ID to 10 (one more than max pattern ID).    */
            constexpr inline int LOOP_ID = 10;

            /* If both deltas in X and Y for a line of specified slope is less than */
            /* this threshold, then the angle for the line is set to 0 radians.     */
            constexpr inline double MIN_SLOPE_DELTA = 0.5;
        }
    }
}
