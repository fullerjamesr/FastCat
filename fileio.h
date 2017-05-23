#pragma once
#ifndef FASTCAT_FILEIO_H
#define FASTCAT_FILEIO_H

#include <stdio.h>
#include <stdbool.h>

size_t readline_resizable(char **buffer, size_t *size, FILE *source, bool incl_newline);
size_t readline(char* buffer, size_t buffer_leng, FILE* source, bool incl_newline);
size_t read_3_column_data(double** buffer, size_t nrows, FILE* source);

#endif //FASTCAT_FILEIO_H