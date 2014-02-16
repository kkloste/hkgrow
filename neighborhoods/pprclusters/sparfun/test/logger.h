/** 
 * @file logger.h
 * Prototypes for the sparvec logger routines
 */

/*
 * David F. Gleich
 * Copyright 2010
 */
 
class Logger
{
public:
    void push_level(const char* name);
    void pop_level();
    
    void log_double(const char* variable, double value);
    void log_string(const char* variable, const char* value);
    void log_int(const char* variable, int value);
    void log_message(const char* format, ...);
    
    
    
    
};

/* usage
 * Logger(
