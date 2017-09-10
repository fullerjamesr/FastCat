#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "argparse.h"


void print_usage(const ArgParseInfo* argParseInfo)
{
    // Print header
    printf("\nUsage: %s\n\n", argParseInfo->usage ? argParseInfo->usage : "");
    if(argParseInfo->description)
        printf("%s\n\n", argParseInfo->description);
    puts("Mandatory arguments to long options are mandatory for short options too.");
    
    // Determine the width of the longest long option
    size_t long_option_len = 0;
    for(const ArgumentOption* option = argParseInfo->options; option && option->type != ARGPARSE_OPT_END; option++)
    {
        if(option->long_name)
        {
            size_t this_len = strlen(option->long_name) + strlen(ArgumentDescriptionStrings[option->type]) + 1;
            if(this_len > long_option_len)
                long_option_len = this_len;
        }
    }
    
    // Describe options
    for(const ArgumentOption* option = argParseInfo->options; option && option->type != ARGPARSE_OPT_END; option++)
    {
        if(option->short_name)
            printf("  -%c", option->short_name);
        else
            fputs("    ", stdout);
        
        if(option->short_name && option->long_name)
            fputs(", ", stdout);
        else
            fputs("  ", stdout);
        
        if(option->long_name)
        {
            char buffer[128] = "";
            strcat(buffer, option->long_name);
            if(strlen(ArgumentDescriptionStrings[option->type]) > 0)
            {
                strcat(buffer, "=");
                strcat(buffer, ArgumentDescriptionStrings[option->type]);
            }
            printf("--%-*s", (int)long_option_len, buffer);
        }
        else
        {
            printf("%-*s  ", (int)long_option_len, ArgumentDescriptionStrings[option->type]);
        }
        
        printf("  %s\n", option->desc ? option->desc : "");
    }
    fputs("\n", stdout);
    
    if(argParseInfo->epilog)
        printf("%s\n", argParseInfo->epilog);
}

void print_state(const ArgParseInfo* argParseInfo)
{
    // Determine the width of the longest option
    size_t long_option_len = 0;
    for(const ArgumentOption* option = argParseInfo->options; option && option->type != ARGPARSE_OPT_END; option++)
    {
        if(!option->value)
            continue;
        if(option->long_name)
        {
            size_t this_len = strlen(option->long_name);
            if(this_len > long_option_len)
                long_option_len = this_len;
        }
    }
    
    for(const ArgumentOption* option = argParseInfo->options; option && option->type != ARGPARSE_OPT_END; option++)
    {
        // Any flag that doesn't have a modifiable value is not the purpose of this function
        if(!option->value)
            continue;
        
        if(option->long_name)
        {
            printf("--%-*s   ", (int)long_option_len, option->long_name);
        }
        else if(option->short_name)
        {
            printf("-%-*c    ", (int)long_option_len, option->short_name);
        }
        
        switch(option->type)
        {
            case ARGPARSE_OPT_BOOLEAN:
            {
                puts(*((bool*)option->value) ? "TRUE" : "FALSE");
                break;
            }
            case ARGPARSE_OPT_STRING:
            {
                puts(*((char**)option->value));
                break;
            }
            case ARGPARSE_OPT_INTEGER:
            {
                printf("%d\n", *((int*)option->value));
                break;
            }
            case ARGPARSE_OPT_SIZE:
            {
                printf("%zu\n", *((size_t*)option->value));
                break;
            }
            case ARGPARSE_OPT_DOUBLE:
            {
                printf("%0.8f\n", *((double*)option->value));
                break;
            }
            default:
            {
                putchar('\n');
            }
        }
    }
}
int do_argparse(const ArgParseInfo* argParseInfo, int argc, const char** argv)
{
    int inputs_count = 0;
    
    for(int i = 1; i < argc; i++)
    {
        const char* token = argv[i];
        
        // Parameter (doesn't start with '-' or is only '-' or is only '--')
        if(token[0] != '-' || !token[1] || (!token[2] && token[1] == '-'))
        {
            argv[inputs_count++] = token;
        }
        // Short option
        else if(token[1] != '-')
        {
            for(size_t j = 1; token[j]; j++)
            {
                const ArgumentOption* option = argParseInfo->options;
                for(; option && option->type != ARGPARSE_OPT_END; option++)
                {
                    if(option->short_name && option->short_name == token[j])
                    {
                        // If this is not a boolean, it can't be combined with anything else
                        if(option->type != ARGPARSE_OPT_BOOLEAN && token[2] != '\0')
                            goto booleancomboerror;
                        
                        char* endptr;
                        switch(option->type)
                        {
                            case ARGPARSE_OPT_HELP:
                            {
                                goto help;
                            }
                            case ARGPARSE_OPT_BOOLEAN:
                            {
                                *( (bool*)option->value ) = *( (bool*)option->value ) ? false : true;
                                break;
                            }
                            case ARGPARSE_OPT_STRING:
                            {
                                if(i + 1 >= argc)
                                    goto requiresvalue;
                                *( (const char**)option->value ) = argv[++i];
                                break;
                            }
                            case ARGPARSE_OPT_INTEGER:
                            {
                                if(i + 1 >= argc)
                                    goto requiresvalue;
                                long int long_value = strtol(argv[++i], &endptr, 0);
                                if(*endptr || long_value > (long int)INT_MAX)
                                    goto parseerror;
                                *( (int*)option->value ) = (int)long_value;
                                break;
                            }
                            case ARGPARSE_OPT_SIZE:
                            {
                                if(i + 1 >= argc)
                                    goto requiresvalue;
                                unsigned long long int ull_value = strtoull(argv[++i], &endptr, 0);
                                if(*endptr || ull_value > (unsigned long long int)SIZE_MAX)
                                    goto parseerror;
                                *( (size_t*)option->value ) = (size_t)ull_value;
                                break;
                            }
                            case ARGPARSE_OPT_DOUBLE:
                            {
                                if(i + 1 >= argc)
                                    goto requiresvalue;
                                double dblvalue = strtod(argv[++i], &endptr);
                                if(*endptr)
                                    goto parseerror;
                                *( (double*)option->value ) = dblvalue;
                                break;
                            }
                            default:
                            {
                                fprintf(stderr, "error: internal argparse configuration error for `%s`\n", token);
                                exit(1);
                            }
                        }
                        break;
                    }
                }
                // This should only be true if we never found a match and thus never broke the above loop early
                if(!option || option->type == ARGPARSE_OPT_END)
                    goto unknown;
            }
        }
        // Long option
        else
        {
            // Find the longest long_name that is a substring of token
            const char* stripped_token = token + 2;
            const ArgumentOption* best_match = NULL;
            for(const ArgumentOption *option = argParseInfo->options; option && option->type != ARGPARSE_OPT_END; option++)
            {
                if(option->long_name && strstr(stripped_token, option->long_name) == stripped_token)
                    if(!best_match || strlen(option->long_name) > strlen(best_match->long_name))
                        best_match = option;
            }
            
            // Did the search find anything at all?
            if(!best_match)
                goto unknown;
            
            if(best_match->type == ARGPARSE_OPT_HELP)
                goto help;
            
            // determine if there's an '=' in this token after the name, or the value is in the next argv
            const char* value = NULL;
            if(strcmp(best_match->long_name, stripped_token) == 0)
            {
                if(best_match->type != ARGPARSE_OPT_BOOLEAN)
                {
                    if((i + 1) < argc)
                        value = argv[++i];
                    else
                        goto requiresvalue;
                }
            }
            else if(best_match->type != ARGPARSE_OPT_BOOLEAN && stripped_token[strlen(best_match->long_name)] == '=')
                value = stripped_token + strlen(best_match->long_name) + 1;
            else
                goto unknown;
            
            char* endptr;
            switch(best_match->type)
            {
                case ARGPARSE_OPT_BOOLEAN:
                {
                    *( (bool*)best_match->value ) = *( (bool*)best_match->value ) ? false : true;
                    break;
                }
                case ARGPARSE_OPT_STRING:
                {
                    *( (const char**) best_match->value) = value;
                    break;
                }
                case ARGPARSE_OPT_INTEGER:
                {
                    long int long_value = strtol(value, &endptr, 0);
                    if(*endptr || long_value > (long int)INT_MAX)
                        goto parseerror;
                    *( (int*) best_match->value ) = (int)long_value;
                    break;
                }
                case ARGPARSE_OPT_SIZE:
                {
                    unsigned long long int ull_value = strtoull(value, &endptr, 0);
                    if(*endptr || ull_value > (unsigned long long int)SIZE_MAX)
                        goto parseerror;
                    *( (size_t*)best_match->value ) = (size_t)ull_value;
                    break;
                }
                case ARGPARSE_OPT_DOUBLE:
                {
                    double dblvalue = strtod(value, &endptr);
                    if(*endptr)
                        goto parseerror;
                    *( (double*) best_match->value ) = dblvalue;
                    break;
                }
                default:
                {
                    fprintf(stderr, "error: internal argparse configuration error for `%s`\n", token);
                    exit(1);
                }
            }
        }
        
        // don't execute through to error conditions
        continue;
        
        // shared fast exit states
        help:
            print_usage(argParseInfo);
            exit(0);
        
        unknown:
            fprintf(stderr, "error: unknown option `%s`\n", token);
            print_usage(argParseInfo);
            exit(1);
        
        booleancomboerror:
            fprintf(stderr, "error: do not combine short flags that require parameters (`%s`)\n", token);
            print_usage(argParseInfo);
            exit(1);
        
        requiresvalue:
            fprintf(stderr, "error: flag `%s` requires value", token);
            print_usage(argParseInfo);
            exit(1);
        
        parseerror:
            fprintf(stderr, "error: could not extract value for `%s` (`%s`?)\n", token, argv[i]);
            print_usage(argParseInfo);
            exit(1);
    }
    return inputs_count;
}
