import sys
import os
import csv
import argparse
import datetime
import logging
import shutil
import yaml

from sarand import params, full_pipeline, utils
from sarand.__init__ import __version__
from sarand.full_pipeline import full_pipeline_main
from sarand.utils import initialize_logger, validate_print_parameters_tools,\
                    update_full_pipeline_params, haveContent

def full_pipeline_init(args, config_data, params):
    """
    """
    # Rewrite parameters of params which are available in data and
    # are supposed to be set for this function
    if config_data=='':
        params = update_full_pipeline_params(args, params)
    else:
        params = update_full_pipeline_params(args, config_data, params)
    #logging file
    log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
    initialize_logger(params.main_dir, log_name)
    #logging.info(str(params.__dict__))

    #create the output directory;
    if not os.path.exists(params.output_dir):
        os.makedirs(params.output_dir)
    #if it exists, delete it and create a new one
    # else:
    #     try:
    #         shutil.rmtree(args.output_dir)
    #     except OSError as e:
    #         logging.error("Error: %s - %s." % (e.filename, e.strerror))
    #     os.makedirs(args.output_dir)
    logging.info('Running full_pipeline ...')
    params = validate_print_parameters_tools(params)
    full_pipeline_main(params)

def main():
    parser = argparse.ArgumentParser(description="Extract the neighborhood of the "
                                        "target Antimicrobial Resistance (AMR) "
                                        "genes from the assembly graph.",
                                        prog='sarand',
                                        usage='sarand <options>')
    parser = create_arguments(parser)
    args = parser.parse_args()
    #If no argument has been passed
    if not len(sys.argv) > 1:
        print("Please use -h option to access usage information!")
        sys.exit()
    # Check if the config file has correctly been set and exist???
    config_data = ''
    if args.config_file!='':
        if not os.path.isfile(args.config_file) or\
            not args.config_file.lower().endswith('yaml') or\
            not haveContent(args.config_file):
            print(args.config_file+" Config file ("+args.config_file+") doesn't exist or not in the right format or is empty! Please provide the correct path to the YAML config file!")
            sys.exit()
        else:
            # Read config file into a dictionery
            print("Reading the config file '"+args.config_file+"' ...")
            with open(args.config_file, 'r') as yamlfile:
                config_data = yaml.load(yamlfile, Loader=yaml.FullLoader)
    #update config info with the ones provided in args
    #data = update_config_info(args, data)

    # if args.config_file=='' or not os.path.isfile(args.config_file) or\
    #     not args.config_file.lower().endswith('yaml'):
    #     print(args.config_file+" doesn't exist or not in the right format! Please provide the correct path to the YAML config file!")
    #     sys.exit()
    # # Read config file into a dictionery
    # print("Reading the config file '"+args.config_file+"' ...")
    # with open(args.config_file, 'r') as yamlfile:
    #     data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    full_pipeline_init(args, config_data, params)
    #args.func(data, params)


if __name__ == '__main__':
    main()
