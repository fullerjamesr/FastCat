#pragma once
#ifndef FASTCAT_ARGPARSE_H
#define FASTCAT_ARGPARSE_H

#include <stddef.h>

typedef enum ArgParseOptionType
{
    ARGPARSE_OPT_END,
    ARGPARSE_OPT_HELP,
    ARGPARSE_OPT_BOOLEAN,
    ARGPARSE_OPT_INTEGER,
    ARGPARSE_OPT_SIZE,
    ARGPARSE_OPT_DOUBLE,
    ARGPARSE_OPT_STRING,
    ARGPARSE__OPT_TYPECOUNT
} ArgParseOptionType;

static const char ArgumentDescriptionStrings[7][15] =
{
        "",
        "",
        "",
        "<int>",
        "<positive int>",
        "<float>",
        "<str>"
};

typedef struct ArgumentOption
{
    ArgParseOptionType type;
    const char short_name;
    const char* long_name;
    void* value;
    const char* desc;
} ArgumentOption;

typedef struct ArgParseInfo
{
    const ArgumentOption* options;
    const char* usage;
    const char* description;
    const char* epilog;
} ArgParseInfo;

// Macros for creating ArgumentOptions conveniently
#define OPT_END()           { ARGPARSE_OPT_END, 0, NULL, NULL, NULL }
#define OPT_HELP(MSG)       { ARGPARSE_OPT_HELP, 'h', "help", NULL, MSG }
#define OPT_BOOLEAN(...)    { ARGPARSE_OPT_BOOLEAN, __VA_ARGS__ }
#define OPT_INTEGER(...)    { ARGPARSE_OPT_INTEGER, __VA_ARGS__ }
#define OPT_SIZE(...)       { ARGPARSE_OPT_SIZE, __VA_ARGS__ }
#define OPT_DOUBLE(...)     { ARGPARSE_OPT_DOUBLE, __VA_ARGS__ }
#define OPT_STRING(...)     { ARGPARSE_OPT_STRING, __VA_ARGS__ }

void print_state(const ArgParseInfo* argParseInfo);
void print_usage(const ArgParseInfo* argParseInfo);
// I don't like the use of `int` here over `size_t` but I use it for consistency with `argc` being an int
int do_argparse(const ArgParseInfo* argParseInfo, int argc, const char** argv);

#endif //FASTCAT_ARGPARSE_H
