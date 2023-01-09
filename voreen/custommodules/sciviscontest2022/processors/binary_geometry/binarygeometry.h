/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#ifndef VRN_BINARYGEOMETRY_H
#define VRN_BINARYGEOMETRY_H

#include <cstdint>
#include <memory>
#include <string>

#ifdef WIN32
#include <winnt.h>
#endif // WIN32

#include "voreen/core/datastructures/geometry/geometry.h"

namespace voreen {

struct GeometryParser;
struct GeometryParserGroup;

class BinaryFileParser {
public:
	static std::unique_ptr<Geometry> read_file(const char* path);
	static void write_file(const Geometry* geometry, const char* path);

	static const GeometryParserGroup* parser_groups() noexcept;
	static std::size_t num_parser_groups() noexcept;

	static const GeometryParser* find_default_parser(const char* type) noexcept;
	static const GeometryParser* find_parser(const char* type, std::uint8_t version) noexcept;
};

class BinaryFile {
	friend BinaryFileParser;
public:
	static constexpr std::size_t FAILED = -1;

	~BinaryFile();

	std::size_t length() const;
	const void* data(std::size_t offset, std::size_t length) const;
	void* data(std::size_t offset, std::size_t length);

	std::size_t align_to(std::size_t offset, std::size_t alignment) const;
	std::size_t align_to(std::size_t offset, std::size_t alignment);

	std::size_t extend_file(std::size_t extra_length);
	void flush();

	template<typename T>
	T* get(std::size_t offset) noexcept {
		return static_cast<T*>(this->data(offset, sizeof(T)));
	}

	template<typename T>
	const T* get(std::size_t offset) const noexcept {
		return static_cast<const T*>(this->data(offset, sizeof(T)));
	}

	template<typename T>
	T* get_array(std::size_t offset, std::size_t count) noexcept {
		return static_cast<T*>(this->data(offset, count * sizeof(T)));
	}

	template<typename T>
	const T* get_array(std::size_t offset, std::size_t count) const noexcept {
		return static_cast<const T*>(this->data(offset, count * sizeof(T)));
	}

	template<typename T>
	std::size_t align_for(std::size_t offset) const noexcept {
		return this->align_to(offset, alignof(T));
	}

	template<typename T>
	std::size_t align_for(std::size_t offset) noexcept {
		return this->align_to(offset, alignof(T));
	}

	template<typename T>
	std::size_t read(T& tmp, std::size_t offset) const noexcept {
		const T* ptr = this->get<T>(offset);
		if (!ptr) {
			return BinaryFile::FAILED;
		}

		std::memcpy(&tmp, ptr, sizeof(T));
		return offset + sizeof(T);
	}

	template<typename T>
	std::size_t read_array(T& tmp, std::size_t offset, std::size_t count) const noexcept {
		const T* ptr = this->get_array<T>(offset, count);
		if (!ptr) {
			return BinaryFile::FAILED;
		}

		std::memcpy(&tmp, ptr, count * sizeof(T));
		return offset + (count * sizeof(T));
	}

	template<typename T>
	std::size_t write(const T& v, std::size_t offset) noexcept {
		T* ptr = this->get<T>(offset);
		if (!ptr) {
			return BinaryFile::FAILED;
		}

		std::memcpy(ptr, &v, sizeof(T));
		return offset + sizeof(T);
	}

	template<typename T>
	std::size_t write_array(const T& v, std::size_t offset, std::size_t count) noexcept {
		T* ptr = this->get_array<T>(offset, count);
		if (!ptr) {
			return BinaryFile::FAILED;
		}

		std::memcpy(ptr, &v, count * sizeof(T));
		return offset + (count * sizeof(T));
	}

	template<typename T>
	std::size_t extend(const T& v) noexcept {
		auto offset = this->align_for<T>(this->length());
		this->extend_file(sizeof(T));
		return this->write(v, offset);
	}

	template<typename T>
	std::size_t extend_array(const T& v, std::size_t count) noexcept {
		auto offset = this->align_for<T>(this->length());
		this->extend_file(count * sizeof(T));
		return this->write_array(v, offset, count);
	}
private:
	BinaryFile() noexcept = default;

	static std::unique_ptr<const BinaryFile> open_read(const char* path);
	static std::unique_ptr<BinaryFile> open_write(const char* path);

	std::size_t file_length;
#ifdef WIN32
	LPVOID mapped_addr;
	HANDLE file;
	HANDLE mapping;
	DWORD protect;
	DWORD access;
#else
	void* mapped_addr;
	int fd;
	int prot;
	int flags;
#endif // WIN32

	const std::string loggerCat_;
};

struct GeometryParserGroup {
	std::string geometry_type;
	const GeometryParser* parsers;
	const GeometryParser* default_parser;
	std::size_t num_parsers;
};

struct GeometryParser {
	std::uint8_t version;
	std::size_t data_alignment;
	std::unique_ptr<Geometry>(*read_file)(const BinaryFile& file, std::size_t offset);
	void(*write_file)(const Geometry* geometry, BinaryFile& file, std::size_t offset);
};

} // namespace

#endif // VRN_BINARYGEOMETRY_H
