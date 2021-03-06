/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.4
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE PACKAGE
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#ifdef PACKAGE_NAME
#define CMDLINE_PARSER_PACKAGE_NAME PACKAGE_NAME
#else
#define CMDLINE_PARSER_PACKAGE_NAME PACKAGE
#endif
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION VERSION
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int int_opt_arg;	/**< @brief Anzahl der Messwerte in den Kanälen.  */
  char * int_opt_orig;	/**< @brief Anzahl der Messwerte in den Kanälen original value given at command line.  */
  const char *int_opt_help; /**< @brief Anzahl der Messwerte in den Kanälen help description.  */
  int frq_opt_arg;	/**< @brief Binäre Kodierung der Indize darzustellender Amplitudenspektren.  */
  char * frq_opt_orig;	/**< @brief Binäre Kodierung der Indize darzustellender Amplitudenspektren original value given at command line.  */
  const char *frq_opt_help; /**< @brief Binäre Kodierung der Indize darzustellender Amplitudenspektren help description.  */
  int time_opt_arg;	/**< @brief Binäre Kodierung der Indize darzustellender Parameter.  */
  char * time_opt_orig;	/**< @brief Binäre Kodierung der Indize darzustellender Parameter original value given at command line.  */
  const char *time_opt_help; /**< @brief Binäre Kodierung der Indize darzustellender Parameter help description.  */
  int samplingrate_opt_arg;	/**< @brief Ganzzahlige Abtastrate in Hz.  */
  char * samplingrate_opt_orig;	/**< @brief Ganzzahlige Abtastrate in Hz original value given at command line.  */
  const char *samplingrate_opt_help; /**< @brief Ganzzahlige Abtastrate in Hz help description.  */
  int channel_start_opt_arg;	/**< @brief Startposition der Kanäle in den Messwerten.  */
  char * channel_start_opt_orig;	/**< @brief Startposition der Kanäle in den Messwerten original value given at command line.  */
  const char *channel_start_opt_help; /**< @brief Startposition der Kanäle in den Messwerten help description.  */
  int channel_count_opt_arg;	/**< @brief Anzahl der Kanäle.  */
  char * channel_count_opt_orig;	/**< @brief Anzahl der Kanäle original value given at command line.  */
  const char *channel_count_opt_help; /**< @brief Anzahl der Kanäle help description.  */
  int qrs_start_opt_arg;	/**< @brief Startposition eines QRS-Komplexes.  */
  char * qrs_start_opt_orig;	/**< @brief Startposition eines QRS-Komplexes original value given at command line.  */
  const char *qrs_start_opt_help; /**< @brief Startposition eines QRS-Komplexes help description.  */
  int qrs_length_opt_arg;	/**< @brief Länge des QRS-Komplexes.  */
  char * qrs_length_opt_orig;	/**< @brief Länge des QRS-Komplexes original value given at command line.  */
  const char *qrs_length_opt_help; /**< @brief Länge des QRS-Komplexes help description.  */
  char * enum_opt_arg;	/**< @brief Eine Zeichenkettenoption mit einer Liste von Werten (default='rms').  */
  char * enum_opt_orig;	/**< @brief Eine Zeichenkettenoption mit einer Liste von Werten original value given at command line.  */
  const char *enum_opt_help; /**< @brief Eine Zeichenkettenoption mit einer Liste von Werten help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int int_opt_given ;	/**< @brief Whether int-opt was given.  */
  unsigned int frq_opt_given ;	/**< @brief Whether frq-opt was given.  */
  unsigned int time_opt_given ;	/**< @brief Whether time-opt was given.  */
  unsigned int samplingrate_opt_given ;	/**< @brief Whether samplingrate-opt was given.  */
  unsigned int channel_start_opt_given ;	/**< @brief Whether channel_start-opt was given.  */
  unsigned int channel_count_opt_given ;	/**< @brief Whether channel_count-opt was given.  */
  unsigned int qrs_start_opt_given ;	/**< @brief Whether qrs_start-opt was given.  */
  unsigned int qrs_length_opt_given ;	/**< @brief Whether qrs_length-opt was given.  */
  unsigned int enum_opt_given ;	/**< @brief Whether enum-opt was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);

extern const char *cmdline_parser_enum_opt_values[];  /**< @brief Possible values for enum-opt. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
