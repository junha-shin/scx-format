#define _GNU_SOURCE
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
static void scalar_path(SEXP x) {
    if (TYPEOF(x) != STRSXP || XLENGTH(x) != 1 || STRING_ELT(x, 0) == NA_STRING)
        Rf_error("path must be one nonmissing string");
}
#ifdef _WIN32
#include <windows.h>
#include <wchar.h>
static wchar_t *wide_path(SEXP x) {
    const char *s = Rf_translateCharUTF8(STRING_ELT(x, 0));
    int n = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, s, -1, NULL, 0);
    if (!n) Rf_error("invalid UTF-8 path");
    wchar_t *w = (wchar_t *) R_alloc(n, sizeof(wchar_t));
    MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, s, -1, w, n);
    return w;
}
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#ifdef __linux__
#include <sys/syscall.h>
#endif
#ifdef __APPLE__
#include <sys/stdio.h>
#endif
#endif

static int is_link_path(SEXP path) {
#ifdef _WIN32
    DWORD attributes = GetFileAttributesW(wide_path(path));
    return attributes != INVALID_FILE_ATTRIBUTES && (attributes & FILE_ATTRIBUTE_REPARSE_POINT);
#else
    struct stat info;
    return lstat(Rf_translateCharUTF8(STRING_ELT(path, 0)), &info) == 0 && S_ISLNK(info.st_mode);
#endif
}
static SEXP scx_is_link(SEXP path) {
    scalar_path(path);
    return Rf_ScalarLogical(is_link_path(path));
}
static SEXP scx_unlink_link(SEXP path) {
    scalar_path(path);
    if (!is_link_path(path)) Rf_error("refusing to remove a non-link path");
#ifdef _WIN32
    wchar_t *wide = wide_path(path);
    DWORD attributes = GetFileAttributesW(wide);
    BOOL success = attributes & FILE_ATTRIBUTE_DIRECTORY ? RemoveDirectoryW(wide) : DeleteFileW(wide);
    if (!success) Rf_error("could not remove link (Windows error %lu)", (unsigned long)GetLastError());
#else
    if (unlink(Rf_translateCharUTF8(STRING_ELT(path, 0))) != 0) Rf_error("could not remove link: %s", strerror(errno));
#endif
    return R_NilValue;
}
static SEXP scx_rename_new(SEXP source, SEXP target) {
    scalar_path(source); scalar_path(target);
#ifdef _WIN32
    if (!MoveFileExW(wide_path(source), wide_path(target), 0))
        Rf_error("no-replace rename failed (Windows error %lu)", (unsigned long) GetLastError());
#elif defined(__linux__) && defined(SYS_renameat2)
    if (syscall(SYS_renameat2, AT_FDCWD, Rf_translateCharUTF8(STRING_ELT(source, 0)),
                AT_FDCWD, Rf_translateCharUTF8(STRING_ELT(target, 0)), 1) != 0)
        Rf_error("no-replace rename failed: %s", strerror(errno));
#elif defined(__APPLE__)
    if (renamex_np(Rf_translateCharUTF8(STRING_ELT(source, 0)),
                   Rf_translateCharUTF8(STRING_ELT(target, 0)), RENAME_EXCL) != 0)
        Rf_error("no-replace rename failed: %s", strerror(errno));
#else
    Rf_error("safe no-replace directory rename unsupported on this platform");
#endif
    return R_NilValue;
}

typedef struct { FILE *file; z_stream stream; int initialized; } gzip_context;
static void gzip_cleanup(void *data) {
    gzip_context *ctx = data;
    if (ctx->initialized) inflateEnd(&ctx->stream);
    if (ctx->file) fclose(ctx->file);
}
static SEXP gzip_body(void *data) {
    gzip_context *ctx = data;
    unsigned char input[65536], output[65536];
    int between = 0, members = 0;
    if (inflateInit2(&ctx->stream, 31) != Z_OK) Rf_error("cannot initialize gzip reader");
    ctx->initialized = 1;
    for (;;) {
        R_CheckUserInterrupt();
        if (ctx->stream.avail_in == 0) {
            size_t n = fread(input, 1, sizeof(input), ctx->file);
            if (!n) {
                if (ferror(ctx->file)) Rf_error("gzip read failed");
                if (!between || !members) Rf_error("truncated gzip stream or missing trailer");
                return R_NilValue;
            }
            ctx->stream.next_in = input;
            ctx->stream.avail_in = (uInt)n;
        }
        if (between) {
            /* Match gzip's concatenated members and zero padding convention. */
            while (ctx->stream.avail_in && *ctx->stream.next_in == 0) {
                ctx->stream.next_in++; ctx->stream.avail_in--;
            }
            if (!ctx->stream.avail_in) continue;
            Bytef *next = ctx->stream.next_in; uInt available = ctx->stream.avail_in;
            if (inflateReset2(&ctx->stream, 31) != Z_OK) Rf_error("cannot reset gzip reader");
            ctx->stream.next_in = next; ctx->stream.avail_in = available;
            between = 0;
        }
        ctx->stream.next_out = output;
        ctx->stream.avail_out = sizeof(output);
        int status = inflate(&ctx->stream, Z_NO_FLUSH);
        if (status == Z_STREAM_END) { between = 1; members++; }
        else if (status != Z_OK && status != Z_BUF_ERROR)
            Rf_error("invalid gzip stream/CRC: %s", ctx->stream.msg ? ctx->stream.msg : "decompression error");
    }
}
static SEXP scx_check_gzip(SEXP path) {
    scalar_path(path);
    gzip_context ctx;
    memset(&ctx, 0, sizeof(ctx));
#ifdef _WIN32
    ctx.file = _wfopen(wide_path(path), L"rb");
#else
    ctx.file = fopen(Rf_translateCharUTF8(STRING_ELT(path, 0)), "rb");
#endif
    if (!ctx.file) Rf_error("cannot open gzip file: %s", strerror(errno));
    return R_ExecWithCleanup(gzip_body, &ctx, gzip_cleanup, &ctx);
}
static const R_CallMethodDef calls[] = {
    {"scx_is_link", (DL_FUNC) &scx_is_link, 1},
    {"scx_unlink_link", (DL_FUNC) &scx_unlink_link, 1},
    {"scx_rename_new", (DL_FUNC) &scx_rename_new, 2},
    {"scx_check_gzip", (DL_FUNC) &scx_check_gzip, 1},
    {NULL, NULL, 0}
};
void attribute_visible R_init_scxSeurat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, calls, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
