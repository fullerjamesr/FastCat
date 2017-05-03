#pragma once
#ifndef FASTCAT_ARGPARSE_H
#define FASTCAT_ARGPARSE_H

#include <stddef.h>

typedef enum ArgumentType
{
    ARGPARSE_OPT_END,
    ARGPARSE_OPT_HELP,
    ARGPARSE_OPT_BOOLEAN,
    ARGPARSE_OPT_INTEGER,
    ARGPARSE_OPT_DOUBLE,
    ARGPARSE_OPT_STRING,
} ArgumentType;

static const char ArgumentDescriptionStrings[6][8] =
{
        "",
        "",
        "",
        "<int>",
        "<float>",
        "<str>"
};

typedef struct ArgumentOption
{
    ArgumentType type;
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
#define OPT_DOUBLE(...)     { ARGPARSE_OPT_DOUBLE, __VA_ARGS__ }
#define OPT_STRING(...)     { ARGPARSE_OPT_STRING, __VA_ARGS__ }


void print_state(ArgParseInfo* argParseInfo);
void print_usage(ArgParseInfo* argParseInfo);
// I don't like the use of `int` here over `size_t` but I use it for consistency with `argc` being an int
int do_argparse(ArgParseInfo* argParseInfo, int argc, const char** argv);

#endif //FASTCAT_ARGPARSE_H
