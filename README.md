# OpenLQM <img src="OpenLQM/openlqm_cmake_utils/nist_spo_two_color.svg" align="right" alt="NIST Information Technology Laboratory" style="width:250px;" />

About
-----
Open Latent Quality Metric, or "OpenLQM," analyzes friction ridge images and returns latent quality metrics.

History
-------
OpenLQM was developed by Noblis, Inc. under contract to the National
Institute of Standards and Technology in 2024-2025, as a port of "LQMetric,"
which was developed by Noblis, Inc. under contract to the Federal Bureau of
Investigation's Criminal Justice Information Services Division in 2012-2014.
More information about LQMetric can be found in the article:

Kalka N. D., Beachler M., Hicklin R. A. "[LQMetric: A Latent Fingerprint Quality Metric for Predicting AFIS Performance and Assessing the Value of Latent Fingerprints]," Journal of Forensic Identification, vol. 70, no. 4, pp. 443–463, Oct 2020.

Download
--------
Pre-built versions of the OpenLQM standalone executable and library for many
platforms are available to download on the
[GitHub Releases](https://github.com/usnistgov/openlqm/releases) page.

Build Dependencies
------------------

Building OpenLQM requires the following dependencies:

 * [OpenCV](https://github.com/opencv/opencv) ([Apache 2 License](https://github.com/opencv/opencv/blob/master/LICENSE))

We use vcpkg to gather and build these packages, including their dependencies.

Quick Build
-----------

> [!IMPORTANT]
> Unless you are *actively developing* code for OpenLQM, we **highly** suggest
> you download from [Releases](https://github.com/usnistgov/openlqm/releases)
> instead of attempting to compile.

> [!NOTE]
> Due to vcpkg use for finding dependencies, Internet access is required during
> configuration in the default configuration.

```bash
git clone https://github.com/usnistgov/openlqm.git

# Set location of vcpkg
VCPKG_ROOT=path/to/vcpkg

# Configure build
cd openlqm
mkdir build
cd build

# Build and package
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" # Windows, add: -DCMAKE_CONFIGURATION_TYPE=Release
cmake --build . # Windows, add --config Release
cpack
```

Communication
-------------
Please note that there are **no plans** for NIST to further enhance or change
OpenLQM. It is meant to be an open source snapshot of LQMetric. Development of
latent image quality algorithms will happen at the international standards level
in the development of ISO/IEC 29794-12.

If you found a bug and can provide steps to reliably reproduce it, or if you
have a feature request, please
[open an issue](https://github.com/usnistgov/openlqm/issues). Other
questions may be addressed to the
[project maintainers](mailto:openlqm@list.nist.gov).


License
-------
OpenLQM is released under the Apache 2.0 license. See the
[LICENSE](https://github.com/usnistgov/openlqm/blob/master/LICENSE.txt)
for details.

[LQMetric: A Latent Fingerprint Quality Metric for Predicting AFIS Performance and Assessing the Value of Latent Fingerprints]: https://fingerprint.nist.gov/openlqm/JFi-2020-4-443.pdf
