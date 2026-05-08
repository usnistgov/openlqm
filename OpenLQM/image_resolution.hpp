/*
 * Copyright 2025 Noblis, Inc.
 * Licensed under the Apache License, Version 2.0.
 */

#pragma once

#include <cstddef>

namespace OpenLQM::Detail {

/// Extract horizontal pixel density in dots per meter from raw image bytes.
/// Returns 0 when the format is unrecognized or no reliable density metadata exists.
unsigned int ReadDotsPerMeterX(const unsigned char* data, std::size_t size);

} // namespace OpenLQM::Detail
