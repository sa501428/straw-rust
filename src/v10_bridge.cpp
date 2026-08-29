#include "../../straw/C++/straw_v10.h"
#include <cstdlib>
#include <cstring>
#include <exception>
#include <string>

extern "C" {
struct StrawV10Record { int32_t x, y; float value; };
using StrawV10Callback = void (*)(void *, StrawV10Record);
static thread_local std::string error;
const char *straw_v10_error() { return error.c_str(); }
void *straw_v10_open(const char *path) { try { error.clear(); return new straw_v10::File(path); } catch (const std::exception &e) { error = e.what(); return nullptr; } }
void straw_v10_close(void *p) { delete static_cast<straw_v10::File *>(p); }
size_t straw_v10_chromosome_count(void *p) { return static_cast<straw_v10::File *>(p)->chromosomes().size(); }
int straw_v10_chromosome(void *p, size_t i, const char **name, int32_t *index, int64_t *length) {
    try { auto c = static_cast<straw_v10::File *>(p)->chromosomes().at(i); *name = strdup(c.name.c_str()); *index = c.index; *length = c.length; return 1; } catch (const std::exception &e) { error = e.what(); return 0; }
}
void straw_v10_string_free(const char *p) { free(const_cast<char *>(p)); }
const char *straw_v10_genome(void *p) { try { error = static_cast<straw_v10::File *>(p)->genome(); return error.c_str(); } catch (...) { return ""; } }
size_t straw_v10_resolution_count(void *p, int frag) { return static_cast<straw_v10::File *>(p)->resolutions(frag ? "FRAG" : "BP").size(); }
int32_t straw_v10_resolution(void *p, int frag, size_t i) { return static_cast<straw_v10::File *>(p)->resolutions(frag ? "FRAG" : "BP").at(i); }
size_t straw_v10_norm_count(void *p) { return static_cast<straw_v10::File *>(p)->normalizations().size(); }
const char *straw_v10_norm(void *p, size_t i) { try { error = static_cast<straw_v10::File *>(p)->normalizations().at(i); return error.c_str(); } catch (...) { return ""; } }
size_t straw_v10_attribute_count(void *p) { return static_cast<straw_v10::File *>(p)->attributes().size(); }
const char *straw_v10_attribute(void *p, size_t i, int value) { try { auto a = static_cast<straw_v10::File *>(p)->attributes().at(i); error = value ? a.second : a.first; return error.c_str(); } catch (...) { return ""; } }
int straw_v10_stream(void *p, const char *type, const char *norm, const char *a, const char *b, const char *unit, int32_t resolution, void *context, StrawV10Callback callback) {
    try { error.clear(); static_cast<straw_v10::File *>(p)->stream(type, norm, a, b, unit, resolution, [&](const contactRecord &r) { callback(context, {r.binX, r.binY, r.counts}); }); return 1; } catch (const std::exception &e) { error = e.what(); return 0; }
}
int64_t straw_v10_count(void *p, int32_t resolution, int inter) { try { error.clear(); return static_cast<straw_v10::File *>(p)->countRecords(resolution, inter != 0); } catch (const std::exception &e) { error = e.what(); return -1; } }
}
