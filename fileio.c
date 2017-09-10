#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "fileio.h"

/*
 * Read characters from a file `source` into a C-string `buffer` until the first newline or EOF or error is reached.
 * Behavior is almost identical to `fgets`, except the number of characters written into the destination is returned in
 * order to make a subsequent call to `strlen(buffer)` unnecessary.
 *
 * Writes at most `buffer_len - 1` characters into `buffer`, then appends a terminating null after the last character
 * written.
 *
 * Parameters:
 *   char* buffer: An array of char of at least length `buffer_len`
 *   size_t buffer_len: The maximum number of chars to write to `buffer`
 *   FILE* source: A pointer to a valid FILE object that can be read from
 *   bool incl_newline: When a newline is encountered, should that be written into `buffer`?
 *
 * Returns:
 *   The number of characters that were written into `buffer`, excluding the terminating null character, which is always
 *   written unless `source` has EOF or ERROR set when the function called. Under normal circumstances this should be
 *   identical to the result of calling strlen(buffer) immediately after this function returns.
 *
 * Errors:
 *   If this returns 0 after being called with `buffer_len` > 1, three possibilities exist
 *   (0) there was a blank line
 *   (1) the final character from `source` had already been read (and so the first read attempt returned EOF)
 *   (2) the EOF or error flags was set on `source` already before this function was called
 *   In the case of #2, `buffer` remains unchanged. Cases #0 and #1 will return 0 and have placed a terminating null at
 *   `buffer[0]`.
 *
 *   Note that if EOF or ERROR results after one or more characters are successfully read, `buffer` will contain those
 *   successfully read characters (and a terminating null) and their count.
 */
size_t readline(char *buffer, size_t buffer_len, FILE *source, bool incl_newline)
{
    assert(buffer);
    assert(source);
    if(feof(source) || ferror(source) || buffer_len < 1)
        return 0;
    
    size_t chars_written = 0;
    int next_char;
    
    while(!ferror(source) && !feof(source) && chars_written < (buffer_len-1))
    {
        next_char = fgetc(source);
        
        if(next_char == EOF)
            break;
        if(!incl_newline && next_char == '\n')
            break;
        
        buffer[chars_written++] = (char)next_char;
        
        if(next_char == '\n')
            break;
    }
    buffer[chars_written] = '\0';
    
    return chars_written;
}

/*
 * Read from a file `source` into the C-string pointed to by `buffer` until the first newline or EOF or error is
 * reached. Behavior is similar to `fgets`, except the destination buffer is realloc'd to try to contain the entire line
 * if possible, and the number of characters written into the destination is returned in order to make a subsequent call
 * to `strlen(buffer)` unnecessary.
 *
 * Parameters:
 *   char** buffer: Characters will be copied into the C-string pointed at by `buffer`. This should be allocated on the
 *                  heap so it can be grown as necessary.
 *   size_t* size: A pointer to the initial/current size of the C-string pointed at by `buffer`. This will be updated if
 *                 the buffer is grown.
 *   FILE* source: A pointer to a valid FILE object to read characters from.
 *   bool incl_newline: When a newline is encountered, should that be written into the buffer?
 *
 * Returns:
 *   The number of characters that were written into `buffer`, excluding the terminating null character, which is always
 *   written unless `source` has EOF or ERROR set when the function called. Under normal circumstances this should be
 *   identical to the result of calling strlen(buffer) immediately after this function returns.
 *
 * Errors:
 *   See `readline` above.
 *   If the buffer needs to be grown but a call to `realloc` fails, the function returns whatever was able to be
 *   successfully read up to that point, just as if EOF had been encountered.
 */
size_t readline_resizable(char **buffer, size_t *size, FILE *source, bool incl_newline)
{
    // sanity checks
    assert(source);
    assert(buffer);
    assert((*buffer));
    assert(size);
    size_t dest_size = *size;
    char* str = *buffer;
    if(feof(source) || ferror(source) || dest_size < 1)
        return 0;
    
    size_t chars_written = 0;
    int next_char;
    
    while(!ferror(source) && !feof(source))
    {
        next_char = fgetc(source);
    
        if(next_char == EOF)
            break;
        if(!incl_newline && next_char == '\n')
            break;
        
        // check if the buffer needs to be resized, and grow geometrically
        // also check for the corner case of dest_size==1 so size_t doesn't underflow
        if(dest_size == 1 || (dest_size-2) <= chars_written)
        {
            size_t new_size = (dest_size + 1) * 2;
            char* new_str = realloc(str, new_size);
            if(new_str == NULL)
                break;
            dest_size = new_size;
            str = new_str;
        }
        
        str[chars_written++] = (char)next_char;
        
        if(next_char == '\n')
            break;
    }
    
    str[chars_written] = '\0';
    *buffer = str;
    *size = dest_size;
    return chars_written;
}

/*
 *
 * Read a 3 column data file, either whitespace or comma delimited, of floating point numbers Each column is read into
 * an array. Any line that does not consist of 1 or more floating point values separated by whitespace or commas is
 * ignored. Any missing values are filled with 0.0
 *
 * Parameters:
 *  double** buffer: a two dimensional array of dimension greater than or equal to [3][nrows]
 *  size_t nrows: the maximum number of rows to read from `source` and store in `buffer[column]`
 *  FILE* source: a valid file to read from
 *
 *  Returns:
 *      The number of lines (rows) successfully read, <= nrows
 *
 */
size_t read_3_column_data(double** buffer, size_t nrows, FILE* source)
{
    size_t rows = 0;
    double one, two, three;
    while(!ferror(source) && !feof(source) && rows < nrows)
    {
        // clear prior values
        one = 0.0, two = 0.0, three = 0.0;
        // grab up to three double values
        if(fscanf(source, "%lf%*[, \t]%lf%*[, \t]%lf", &one, &two, &three) > 0)
        {
            buffer[0][rows] = one;
            buffer[1][rows] = two;
            buffer[2][rows++] = three;
        }
        // move the stream ahead until a newline appears
        while(!ferror(source) && !feof(source) && fgetc(source) != '\n');
    }
    return rows;
}

/*
 *
 * Write double values as a tab delimited file. Each array becomes a column.
 *
 * Parameters:
 *  FILE* destination: A valid pointer to a writable file or stream
 *  size_t nrows: The number of rows to write. This should be <= the length of all the data arrays.
 *  size_t ncols: The number of arrays (of type double) being passed, and thus the number of columns to write
 *  ...: Any number (>= `ncols`) of arrays of type double to be written. None should be shorter than `nrows`.
 *
 * Returns:
 *  The number of rows successfully written. Always equal to `nrows` unless an error flag was raised on the stream.
 *
 * Errors:
 *  If the error flag is set on `destination`, either beforehand or in the process of writing, this function will return
 *  the number of rows successfully written. This can be compared to the input `nrows` to see if something went wrong.
 */
size_t write_n_column_data(FILE* destination, const size_t nrows, const size_t ncols,  ...)
{
    assert(destination);
    
    va_list vargs;
    size_t r;
    for(r = 0; r < nrows && !ferror(destination); r++)
    {
        va_start(vargs, ncols);
        for(size_t c = 0; c < ncols && !ferror(destination); c++)
        {
            double* col = va_arg(vargs, double*);
            char followup = (char) (c == ncols - 1 ? '\n' : '\t');
            if(col)
                fprintf(destination, "%lf%c", col[r], followup);
            else
                fputc(followup, destination);
        }
        va_end(vargs);
    }
    return r;
}
