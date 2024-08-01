/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "binarygeometry.h"

#include <limits>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif // WIN32

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

constexpr const char* BINARY_FILE_LOGGER_CAT = "voreen.sciviscontest2022.BinaryFile";
constexpr const char* BINARY_FILE_PARSER_LOGGER_CAT = "voreen.sciviscontest2022.BinaryFileParser";

std::unique_ptr<Geometry> read_pointsegmentlistgeometryvec3(const BinaryFile& file, std::size_t offset);
void write_pointsegmentlistgeometryvec3(const Geometry* geometry, BinaryFile& file, std::size_t offset);

///< Header for a PointSegmentListGeometryVec3 data segment.
struct PointSegmentListGeometryHeader {
	std::uint64_t data_start;
	std::uint64_t num_segments;
	std::uint8_t has_transform_matrix;
};

const GeometryParser point_segment_list_parsers[] = {
	{
		1,
		alignof(PointSegmentListGeometryHeader),
		read_pointsegmentlistgeometryvec3,
		write_pointsegmentlistgeometryvec3
	}
};

const GeometryParserGroup geometry_parser_groups[] = {
	{
		"PointSegmentListGeometryVec3",
		point_segment_list_parsers,
		&point_segment_list_parsers[0],
		sizeof(point_segment_list_parsers)/sizeof(GeometryParser)
	}
};

constexpr std::uint32_t MAGIC_NUMBER = 0x45474256;	///< "VBGE"

struct BinaryGeometryHeader {
	std::uint32_t magic;	///< Magic number.
	std::uint8_t version;	///< Version of the binary.
	std::uint8_t geometry_type_length;	///< Length of the geometry type string.
	std::uint16_t data_segment_offset;	///< Offset from the beginning of the file to the start of the data segment.
	std::uint8_t parser_version;	///< Version of the parser used to create the file.
};

static_assert(sizeof(BinaryGeometryHeader) == 12, "BinaryGeometryHeader must be 12 bytes long");

std::unique_ptr<Geometry> BinaryFileParser::read_file(const char *path) {
	auto file = BinaryFile::open_read(path);
	if (!file) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "file could not be opened");
		return {};
	}

	BinaryGeometryHeader header {};
	auto geometry_type_offset = file->read(header, 0);
	if (geometry_type_offset == BinaryFile::FAILED || header.magic != MAGIC_NUMBER) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "file does not contain a valid geometry header");
		return {};
	}

	// We only define the header for the version 1.
	//
	// File structure
	// 0..4								MAGIC_STRING
	// 4..5 							HEADER_VERSION
	// 5..6 							GEOMETRY_TYPE_LENGTH
	// 6..8 							DATA_SEGMENT_OFFSET
	// 8..9								PARSER_VERSION
	// 9..12							PADDING (Reserved)
	// 12..12 + GEOMETRY_TYPE_LENGTH	GEOMETRY_TYPE
	// ...								PADDING
	// DATA_SEGMENT_OFFSET..END			DATA_SEGMENT
	if (header.version != 1) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "unrecognized header version");
		return {};
	}

	// Fetch geometry type.
	// Geometry type is stored after the header.
	char geometry_type[257] = { '\0' };		// extra byte for '\0' terminator.
	auto geometry_type_end = file->read_array(geometry_type[0], geometry_type_offset, header.geometry_type_length);
	if (geometry_type_end == BinaryFile::FAILED) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "file does not contain a valid geometry type");
		return {};
	}

	if (geometry_type_end > header.data_segment_offset) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "data segment points into the file header");
		return {};
	}

	// find a suitable parser.
	auto parser = BinaryFileParser::find_parser(geometry_type, header.parser_version);
	if (parser == nullptr) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "no suitable parser found, type: " << geometry_type);
		return {};
	}

	// check alignment.
	if (header.data_segment_offset % parser->data_alignment != 0) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "data segment is not aligned");
		return {};
	}

	return parser->read_file(*file, header.data_segment_offset);
}

void BinaryFileParser::write_file(const Geometry *geometry, const char *path) {
	if (geometry==nullptr) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "geometry is null");
		return;
	}

	auto file = BinaryFile::open_write(path);
	if (!file) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "could not open file");
		return;
	}

	auto geometry_type = geometry->getClassName();
	if (geometry_type.length() > 256) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "invalid geometry type length");
		return;
	}

	auto header_ptr = file->get<BinaryGeometryHeader>(0);
	if (!header_ptr || file->length() != sizeof(BinaryGeometryHeader)) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "invalid file");
		return;
	}

	header_ptr->magic = MAGIC_NUMBER;
	header_ptr->version = 1;
	header_ptr->geometry_type_length = static_cast<std::uint8_t>(geometry_type.length());
	auto geometry_type_end = file->extend_array(geometry_type[0], header_ptr->geometry_type_length);

	if (geometry_type_end == BinaryFile::FAILED) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "could not write file header");
		return;
	}

	// find a suitable parser.
	auto parser = BinaryFileParser::find_default_parser(geometry_type.c_str());
	if (parser == nullptr) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "no suitable parser found");
		return;
	}

	// Extend file up to the start of the data segment.
	auto data_segment_start = file->align_to(geometry_type_end, parser->data_alignment);
	if (data_segment_start == BinaryFile::FAILED) {
		LERRORC(BINARY_FILE_PARSER_LOGGER_CAT, "could not initialize data segment");
		return;
	}

	// Write back data segment start and parser version.
	header_ptr = file->get<BinaryGeometryHeader>(0);
	header_ptr->data_segment_offset = data_segment_start;
	header_ptr->parser_version = parser->version;

	// Parse the data segment.
	parser->write_file(geometry, *file, data_segment_start);
}

const GeometryParserGroup* BinaryFileParser::parser_groups() noexcept {
	return geometry_parser_groups;
}

std::size_t BinaryFileParser::num_parser_groups() noexcept {
	return sizeof(geometry_parser_groups)/sizeof(GeometryParserGroup);
}

const GeometryParser* BinaryFileParser::find_default_parser(const char* type) noexcept {
	if (type == nullptr) {
		return nullptr;
	}

	for (const auto& group: geometry_parser_groups) {
		if (group.geometry_type == type) {
			return group.default_parser;
		}
	}

	return nullptr;
}

const GeometryParser* BinaryFileParser::find_parser(const char* type, std::uint8_t version) noexcept {
	if (type == nullptr) {
		return nullptr;
	}

	for (const auto& group: geometry_parser_groups) {
		if (group.geometry_type == type) {
			for (auto i = 0u; i < group.num_parsers; ++i) {
				if (group.parsers[i].version == version) {
					return &group.parsers[i];
				}
			}
		}
	}

	return nullptr;
}

BinaryFile::~BinaryFile() {
#ifdef WIN32
	if (this->mapped_addr != nullptr) {
		UnmapViewOfFile(this->mapped_addr);
		CloseHandle(this->mapping);
		CloseHandle(this->file);

		this->mapped_addr = nullptr;
		this->mapping = nullptr;
		this->file = nullptr;
	}
#else
	if (this->mapped_addr!=MAP_FAILED) {
		munmap(this->mapped_addr, this->file_length);
		close(this->fd);

		this->mapped_addr = nullptr;
		this->fd = 0;
	}
#endif // WIN32
}

std::size_t BinaryFile::length() const {
	return this->file_length;
}

const void *BinaryFile::data(std::size_t offset, std::size_t length) const {
	if (offset + length > this->file_length) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "offset exceeds file length");
		return nullptr;
	}

#ifdef WIN32
	if (this->mapped_addr == nullptr) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "no address mapped");
		return nullptr;
	}

	const char *ptr = static_cast<const char *>(this->mapped_addr);
	return ptr + offset;
#else
	if (this->mapped_addr==MAP_FAILED) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "no address mapped");
		return nullptr;
	}

	const char *ptr = static_cast<const char *>(this->mapped_addr);
	return ptr + offset;
#endif // WIN32
}

void *BinaryFile::data(std::size_t offset, std::size_t length) {
	if (offset + length > this->file_length) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "offset exceeds file length");
		return nullptr;
	}

#ifdef WIN32
	if (this->mapped_addr == nullptr) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "no address mapped");
		return nullptr;
	}

	char *ptr = static_cast<char *>(this->mapped_addr);
	return ptr + offset;
#else
	if (this->mapped_addr==MAP_FAILED) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "no address mapped");
		return nullptr;
	}

	char *ptr = static_cast<char *>(this->mapped_addr);
	return ptr + offset;
#endif // WIN32
}

std::size_t BinaryFile::align_to(std::size_t offset, std::size_t alignment) const {
	if (offset > this->file_length) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "invalid offset");
		return BinaryFile::FAILED;
	}

	auto aligned_offset = (((offset) + (alignment - 1)) & ~(alignment - 1));
	if (aligned_offset > this->file_length) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "alignment operation out of bounds");
		return BinaryFile::FAILED;
	}

	return aligned_offset;
}

std::size_t BinaryFile::align_to(std::size_t offset, std::size_t alignment) {
	if (offset > this->file_length) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "invalid offset");
		return BinaryFile::FAILED;
	}

	auto aligned_offset = (((offset) + (alignment - 1)) & ~(alignment - 1));
	if (aligned_offset > this->file_length) {
		auto required_bytes = aligned_offset - this->file_length;
		return this->extend_file(required_bytes);
	}

	return aligned_offset;
}

std::size_t BinaryFile::extend_file(std::size_t extra_length) {
	std::size_t new_length = this->file_length + extra_length;

#ifdef WIN32
	// Extend file.
	UnmapViewOfFile(this->mapped_addr);
	CloseHandle(this->mapping);

	this->mapped_addr = nullptr;
	this->mapping = nullptr;

	LARGE_INTEGER li;
	li.QuadPart = new_length;

	DWORD result = SetFilePointer(this->file, li.LowPart, &li.HighPart, FILE_BEGIN);
	if (result == INVALID_SET_FILE_POINTER) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "SetFilePointer failed");
		CloseHandle(this->file);
		this->file = nullptr;
		return BinaryFile::FAILED;
	}

	if (!SetEndOfFile(this->file)) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "SetEndOfFile failed");
		CloseHandle(this->file);
		this->file = nullptr;
		return BinaryFile::FAILED;
	}
    this->file_length = new_length;

	// Go to the start of the file.
	result = SetFilePointer(this->file, 0, nullptr, FILE_BEGIN);
	if (result == INVALID_SET_FILE_POINTER) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "SetFilePointer failed");
		CloseHandle(this->file);
		this->file = nullptr;
		return BinaryFile::FAILED;
	}

	// Map file.
	this->mapping = CreateFileMappingA(file, nullptr, this->protect, 0, 0, nullptr);
	if (this->mapping == nullptr) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "CreateFileMappingA failed");
		CloseHandle(this->file);
		this->file = nullptr;
		return BinaryFile::FAILED;
	}

	this->mapped_addr = MapViewOfFile(this->mapping, this->access, 0, 0, 0);
	if (this->mapped_addr == nullptr) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "MapViewOfFile failed");
		CloseHandle(this->mapping);
		CloseHandle(this->file);
		this->mapping = nullptr;
		this->file = nullptr;
		return BinaryFile::FAILED;
	}
#else
	constexpr off_t MAX_LENGTH = std::numeric_limits<off_t>::max();
	if (new_length > MAX_LENGTH) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "maximum file length reached");
		return BinaryFile::FAILED;
	}

	// Extend file.
	if (ftruncate(this->fd, static_cast<off_t>(new_length))==-1) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "ftruncate failed");
		return BinaryFile::FAILED;
	}

	// map the file again.
	if (munmap(this->mapped_addr, this->file_length)==-1) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "munmap failed");
		return BinaryFile::FAILED;
	}

	this->file_length = new_length;
	this->mapped_addr = mmap(nullptr, new_length, this->prot, this->flags, this->fd, 0);
	if (this->mapped_addr==MAP_FAILED) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "mmap failed");
		close(this->fd);
		this->fd = 0;
		return BinaryFile::FAILED;
	}
#endif // WIN32

	return this->file_length;
}

void BinaryFile::flush() {
#ifdef WIN32
	if (this->mapped_addr == nullptr) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "no address mapped");
		return;
	}

	if (!FlushViewOfFile(this->mapped_addr, this->file_length)) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "FlushViewOfFile failed");
	}
#else
	if (this->mapped_addr==MAP_FAILED) {
		LWARNINGC(BINARY_FILE_LOGGER_CAT, "no address mapped");
		return;
	}

	if (msync(this->mapped_addr, this->file_length, MS_SYNC)==-1) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "msync failed");
	}
#endif // WIN32
}

std::unique_ptr<const BinaryFile> BinaryFile::open_read(const char *path) {
#ifdef WIN32
	HANDLE file = CreateFileA(path, GENERIC_READ, 0, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (file == INVALID_HANDLE_VALUE) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "CreateFileA failed");
		return {};
	}

	// Fetch file length.
	ULARGE_INTEGER li {};
	li.LowPart = GetFileSize(file, &li.HighPart);
	if (li.LowPart == INVALID_FILE_SIZE) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "GetFileSize failed");
		CloseHandle(file);
		return {};
	}

	// Map file.
	HANDLE mapping = CreateFileMappingA(file, nullptr, PAGE_READONLY, 0, 0, nullptr);
	if (mapping == nullptr) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "CreateFileMappingA failed");
		CloseHandle(file);
		return {};
	}

	LPVOID addr = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
	if (addr == nullptr) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "MapViewOfFile failed");
		CloseHandle(mapping);
		CloseHandle(file);
		return {};
	}

	std::unique_ptr<BinaryFile> b_file{new BinaryFile()};
    b_file->file = file;
    b_file->file_length = li.QuadPart;
    b_file->mapping = mapping;
    b_file->protect = PAGE_READONLY;
    b_file->access = FILE_MAP_READ;
    b_file->mapped_addr = addr;

	return b_file;
#else
	int fd = open(path, O_RDONLY);
	if (fd==-1) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "open failed");
		return {};
	}

	// Fetch file length.
	struct stat sb{};
	if (fstat(fd, &sb)==-1) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "fstat failed");
		close(fd);
		return {};
	}

	// Map file.
	void *addr = mmap(nullptr, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
	if (addr==MAP_FAILED) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "mmap failed");
		close(fd);
		return {};
	}

	std::unique_ptr<BinaryFile> file{new BinaryFile()};
	file->fd = fd;
	file->file_length = sb.st_size;
	file->prot = PROT_READ;
	file->flags = MAP_SHARED;
	file->mapped_addr = addr;

	return file;
#endif // WIN32
}

std::unique_ptr<BinaryFile> BinaryFile::open_write(const char *path) {
#ifdef WIN32
	HANDLE file = CreateFileA(path, GENERIC_READ | GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (file == INVALID_HANDLE_VALUE) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "CreateFileA failed");
		return {};
	}

	LARGE_INTEGER li;
	li.QuadPart = sizeof(BinaryGeometryHeader);

	// Extend file.
	DWORD result = SetFilePointer(file, li.LowPart, &li.HighPart, FILE_BEGIN);
	if (result == INVALID_SET_FILE_POINTER) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "SetFilePointer failed");
		CloseHandle(file);
		return {};
	}

	if (!SetEndOfFile(file)) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "SetEndOfFile failed");
		CloseHandle(file);
		return {};
	}

	// Go to the start of the file.
	result = SetFilePointer(file, 0, nullptr, FILE_BEGIN);
	if (result == INVALID_SET_FILE_POINTER) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "SetFilePointer failed");
		CloseHandle(file);
		return {};
	}

	// Map file.
	HANDLE mapping = CreateFileMappingA(file, nullptr, PAGE_READWRITE, 0, 0, nullptr);
	if (mapping == nullptr) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "CreateFileMappingA failed");
		CloseHandle(file);
		return {};
	}

	LPVOID addr = MapViewOfFile(mapping, FILE_MAP_READ | FILE_MAP_WRITE, 0, 0, 0);
	if (addr == nullptr) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "MapViewOfFile failed");
		CloseHandle(mapping);
		CloseHandle(file);
		return {};
	}

	std::unique_ptr<BinaryFile> b_file{new BinaryFile()};
    b_file->file = file;
    b_file->file_length = sizeof(BinaryGeometryHeader);
    b_file->mapping = mapping;
    b_file->protect = PAGE_READWRITE;
    b_file->access = FILE_MAP_READ | FILE_MAP_WRITE;
    b_file->mapped_addr = addr;

	return b_file;
#else
	int fd = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
	if (fd==-1) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "open failed");
		return {};
	}

	// Extend file.
	if (ftruncate(fd, sizeof(BinaryGeometryHeader))==-1) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "ftruncate failed");
		return {};
	}

	// Map file.
	void *addr = mmap(nullptr, sizeof(BinaryGeometryHeader), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (addr==MAP_FAILED) {
		LERRORC(BINARY_FILE_LOGGER_CAT, "mmap failed");
		close(fd);
		return {};
	}

	std::unique_ptr<BinaryFile> file{new BinaryFile()};
	file->fd = fd;
	file->file_length = sizeof(BinaryGeometryHeader);
	file->prot = PROT_READ | PROT_WRITE;
	file->flags = MAP_SHARED;
	file->mapped_addr = addr;

	return file;
#endif // WIN32
}

constexpr const char* POINTSEGMENTLISTVEC3_LOGGER_CAT = "voreen.sciviscontest2022.BinaryFileParser.PointSegmentListVec3";

///< PointSegmentListGeometry are represented as linked lists in binary.
struct PointSegmentListGeometryNode {
	std::uint64_t next;
	std::uint64_t data;
	std::uint64_t data_count;
};

std::unique_ptr<Geometry> read_pointsegmentlistgeometryvec3(const BinaryFile& file, std::size_t offset) {
	std::unique_ptr<PointSegmentListGeometryVec3> geometry{ new PointSegmentListGeometryVec3() };

	PointSegmentListGeometryHeader header {};
	if (file.read(header, offset) == BinaryFile::FAILED) {
		LERRORC(POINTSEGMENTLISTVEC3_LOGGER_CAT, "missing data header");
		return {};
	}

	std::size_t next = header.data_start;
	std::size_t data_end = header.data_start;

	// Read segments.
	if (header.num_segments != 0) {
		std::vector<std::vector<tgt::vec3>> geometry_data = {};
		geometry_data.reserve(header.num_segments);

		while (next != 0) {
			auto node_addr = file.get<PointSegmentListGeometryNode>(next);
			if (!node_addr) {
				LERRORC(POINTSEGMENTLISTVEC3_LOGGER_CAT, "missing segment node");
				return {};
			}

			std::vector<tgt::vec3> segment {};
			segment.resize(node_addr->data_count);
			data_end = file.read_array(segment[0], node_addr->data, node_addr->data_count);
			if (data_end == BinaryFile::FAILED) {
				LERRORC(POINTSEGMENTLISTVEC3_LOGGER_CAT, "missing segment data");
				return {};
			}

			geometry_data.push_back(std::move(segment));

			next = node_addr->next;
		}

		geometry->setData(std::move(geometry_data));
	}

	// Read transformation matrix.
	if (header.has_transform_matrix) {
		// First align to PointSegmentListGeometryNode then to a tgt::mat4.
		auto matrix_start = file.align_for<PointSegmentListGeometryNode>(data_end);
		matrix_start = file.align_for<tgt::mat4>(matrix_start);

		tgt::mat4 m;
		if (file.read(m, matrix_start) == BinaryFile::FAILED) {
			LERRORC(POINTSEGMENTLISTVEC3_LOGGER_CAT, "missing transformation matrix");
			return {};
		}
		geometry->setTransformationMatrix(m);
	}

	return geometry;
}

void write_pointsegmentlistgeometryvec3(const Geometry* geometry, BinaryFile& file, std::size_t offset) {
	auto segment_list = dynamic_cast<const PointSegmentListGeometryVec3*>(geometry);
	if (segment_list == nullptr) {
		LERRORC(POINTSEGMENTLISTVEC3_LOGGER_CAT, "invalid geometry type");
		return;
	}

	auto& geometry_data = segment_list->getData();

	PointSegmentListGeometryHeader header {};
	header.num_segments = geometry_data.size();
	header.has_transform_matrix = true;

	PointSegmentListGeometryNode* previous = nullptr;
	auto segments_start = file.extend(header);
	segments_start = file.align_for<PointSegmentListGeometryNode>(segments_start);

	auto header_ptr = file.get<PointSegmentListGeometryHeader>(offset);
	header_ptr->data_start = segments_start;
	auto next = segments_start;

	// Write segments.
	for (const auto& segment: geometry_data) {
		if (previous != nullptr) {
			previous->next = next;
		}

		auto node_addr_end = file.extend<PointSegmentListGeometryNode>({});
		auto node_addr = file.get<PointSegmentListGeometryNode>(next);

		node_addr->data_count = segment.size();
		node_addr->next = 0;
		node_addr->data = file.align_for<tgt::vec3>(node_addr_end);

		auto data_end = file.extend_array<tgt::vec3>(segment[0], node_addr->data_count);

		auto current = next;
		next = file.align_for<PointSegmentListGeometryNode>(data_end);
		previous = file.get<PointSegmentListGeometryNode>(current);
	}

	// Append transformation matrix.
	auto matrix = segment_list->getTransformationMatrix();
	file.extend(matrix);
}

}
