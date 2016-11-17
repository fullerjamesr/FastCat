#include <assert.h>
#include <stdlib.h>
#include "fileio.h"

/*
Read from a file until the first newline or EOF or error is reached.
Behavior similar to fgets, except the destination buffer can be realloc'd.

Parameters:
    char** buffer: The line will be read into `buffer`
    size_t* size: The current size of `buffer`. This will be grown geometrically if buffer is realloc'd
    FILE* source: A pointer to a valid FILE object that can be read from
    bool incl_newline: When a newline is encountered, should that be written into `buffer`?

Returns:
    The number of characters that were written into *buffer, including the CString null terminator.
    If this returns 0, something went wrong.
*/
size_t readline(char** buffer, size_t* size, FILE* source, bool incl_newline)
{
    size_t dest_size = *size;
    char* str = *buffer;
    
    // sanity checks
    assert(source != NULL);
    assert(buffer != NULL);
    assert((*buffer) != NULL);
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
        
        // check if needs resize and grow geometrically
        if((dest_size-2) <= chars_written)
        {
            size_t new_size = (dest_size+1)*2;
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
    
    str[chars_written++] = '\0';
    *buffer = str;
    *size = dest_size;
    return chars_written;
}
