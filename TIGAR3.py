#!/usr/bin/env python

### Import modules
import argparse, multiprocessing, os, sys
from time import time


# import importlib
# importlib.reload(Job)


class Helper_Job():
	def __init__(self):
		pass
	def add_arguments(self, *args):
		pass

### Functions for setting up job
def help_func():
	return Helper_Job()

def tigar_train_dpr():
	import Model_Train_Pred.DPR_Train as Job
	return Job
	
def tigar_train_en():
	import Model_Train_Pred.Elastic_Net_Train as Job
	return Job

def tigar_ld():
	import Get_LD as Job
	return Job

def tigar_bgwtwas():
	import BGWTWAS_Asso_Study as Job
	return Job

def tigar_twas_indv():
	import Asso_Study_01 as Job
	return Job

def tigar_twas_summ():
	import Asso_Study_02 as Job
	return Job

def tigar_vctwas_indv():
	import VC_TWAS_individual as Job
	return Job

def tigar_vctwas_summ():
	import VC_TWAS_summary  as Job
	return Job

def tigar_predict():
	import Prediction as Job
	return Job


if __name__ == '__main__':
	### time calculation
	start_time = time()

	### set up parser
	init_parser = argparse.ArgumentParser()
	init_parser.add_argument('--TIGAR_dir', type=str, help='tool directory', 
		default=os.path.abspath(os.path.dirname(__file__)))

	### set up command subparsers
	subparsers = init_parser.add_subparsers(dest='command')

	## TRAIN
	parser_train = subparsers.add_parser('train', help='Train GReX prediction models', add_help=False)
	parser_train_model = parser_train.add_subparsers(dest='model', help='MODEL: model used for training (dpr: Bayesion Dirichlet Process Regression Model; en,elasticnet,predixcan: Elastic-Net/PrediXcan model)')
	parser_train_dpr = parser_train_model.add_parser('dpr', add_help=False)
	parser_train_en = parser_train_model.add_parser('en', aliases=['elasticnet','predixcan'], add_help=False)
	parser_train.set_defaults(func=help_func)
	parser_train_dpr.set_defaults(func=tigar_train_dpr)
	parser_train_en.set_defaults(func=tigar_train_en)

	## PRED
	parser_predict = subparsers.add_parser('predict', help='Predict GReX', aliases=['pred'], add_help=False)
	parser_predict.set_defaults(func=tigar_predict)

	## LD
	parser_ld = subparsers.add_parser('ld', help='Generate reference LD covariance files for summary-level TWAS', add_help=False)
	parser_ld.set_defaults(func=tigar_ld)

	## BGW-TWAS
	parser_bgwtwas = subparsers.add_parser('bgwtwas', help='BGW-TWAS', add_help=False)
	parser_bgwtwas.set_defaults(func=tigar_bgwtwas)

	## TWAS
	parser_twas = subparsers.add_parser('twas', help='TWAS', add_help=False)
	parser_twas_assoc = parser_twas.add_subparsers(dest='assoc',
		help='ASSOC: association study to do (1,indv: individual level data; 2,summ: summary level data)')
	parser_twas_indv = parser_twas_assoc.add_parser('indv', aliases=['1'], add_help=False)
	parser_twas_summ = parser_twas_assoc.add_parser('summ', aliases=['2'], add_help=False)
	parser_twas.set_defaults(func=help_func)
	parser_twas_indv.set_defaults(func=tigar_twas_indv)
	parser_twas_summ.set_defaults(func=tigar_twas_summ)

	## VC-TWAS
	parser_vctwas = subparsers.add_parser('vctwas', help='VC-TWAS', add_help=False)
	parser_vctwas_assoc = parser_vctwas.add_subparsers(dest='assoc',
		help='ASSOC: association study to do (1,indv: individual level data; 2,summ: summary level data)')
	parser_vctwas_indv = parser_vctwas_assoc.add_parser('indv', aliases=['1'], add_help=False)
	parser_vctwas_summ = parser_vctwas_assoc.add_parser('summ', aliases=['2'], add_help=False)
	parser_vctwas.set_defaults(func=help_func)
	parser_vctwas_indv.set_defaults(func=tigar_vctwas_indv)
	parser_vctwas_summ.set_defaults(func=tigar_vctwas_summ)

	### Parse known args
	args, unkargs = init_parser.parse_known_args()

	### subparser dict
	if (args.command == 'train'):
		if args.model:
			subcommand = args.model
		else:
			subcommand = 'MODEL'
	elif (args.command in ['twas', 'vctwas']):
		if args.assoc:
			subcommand = args.assoc
		else:
			subcommand = 'ASSOC'
	else:
		subcommand = 'NA'

	### subparser dict
	sp_dict = {
		'train' : {'MODEL': parser_train, 'dpr' : parser_train_dpr, 
			'en' : parser_train_en,
			'elasticnet' : parser_train_en,
			'predixcan' : parser_train_en},
		'predict' : {'NA' : parser_predict},
		'pred' : {'NA' : parser_predict},
		'ld' : {'NA' : parser_ld},
		'bgwtwas' : {'NA' : parser_bgwtwas},
		'twas' : {'ASSOC': parser_twas, 'indv' : parser_twas_indv, 'summ': parser_twas_summ, '1' : parser_twas_indv, '2': parser_twas_summ},
		'vctwas' : {'ASSOC': parser_vctwas, 'indv' : parser_vctwas_indv, 'summ': parser_vctwas_summ, '1' : parser_vctwas_indv, '2': parser_vctwas_summ}}

	subcom_dict = {
		'na' : '',
		'model' : '',
		'assoc': '',
		'dpr' : 'using DPR model',
		'en' : 'using Elastic-Net model',
		'elasticnet' : 'using Elastic-Net model',
		'predixcan' : 'using Elastic-Net model',
		'1' : 'using individual-level data',
		'indv' : 'using individual-level data',
		'2' : 'using summary-level data',
		'summ' : 'using summary-level data'}

	usage_dict = {
		'train' : 'train ' + subcommand,
		'predict' : 'predict',
		'pred' : 'pred',
		'ld' : 'ld',
		'bgwtwas' : 'bgwtwas',
		'twas' : 'twas ' + subcommand,
		'vctwas' : 'vctwas ' + subcommand}

	desc_dict = {
		'train' : 'Train GReX prediction models ' + subcom_dict[subcommand.lower()],
		'predict' : 'predict GReX',
		'pred' : 'predict GReX',
		'ld' : 'generate reference LD covariance files for summary-level TWAS',
		'bgwtwas' : 'conduct TWAS with BGW-TWAS weights',
		'twas' : 'conduct TWAS ' + subcom_dict[subcommand.lower()] ,
		'vctwas' : 'conduct VC-TWAS ' + subcom_dict[subcommand.lower()]}

	### get subparser
	sp = sp_dict[args.command][subcommand]
	usage = 'TIGAR.py ' + usage_dict[args.command]
	desc = desc_dict[args.command]

	parser = argparse.ArgumentParser(parents=[sp], usage=usage, description=desc)
	parser.add_argument('--check_input', action='store_true', 
		help='test intended input; check input args, create directories, output to log but do not run job; all required arguments must be supplied')
	
	### Import TIGARutils
	sys.path.append(args.TIGAR_dir)
	sys.path.append(args.TIGAR_dir + '/Model_Train_Pred')
	sys.path.append(args.TIGAR_dir + '/TWAS')
	sys.path.append(args.TIGAR_dir + '/VC_TWAS')

	import TIGARutils as tg

	### set job module
	Job = args.func()

	### add arguments
	Job.add_arguments(parser)

	### parse remaining arguments
	args.__dict__.update(parser.parse_args(unkargs).__dict__)
	args = Job.setup_arguments(args)

	### Make output, log directories
	os.makedirs(args.out_sub_dir, exist_ok=True)
	os.makedirs(os.path.join(args.out_dir, 'logs'), exist_ok=True)

	### Make temporary directories, if applicable
	if args.__dict__.get('tmp_file_dir'):
		os.makedirs(args.tmp_file_dir, exist_ok=True)

	if args.__dict__.get('tmp_file_cv_dir'):
		os.makedirs(args.tmp_file_cv_dir, exist_ok=True)

	### set stdout to log
	sys.stdout = open(os.path.join(args.out_dir, 'logs', args.log_file), 'w')

	### Check tabix command
	tg.check_tabix()

	### Check input files
	tg.check_input_files(args)

	### print arguments
	Job.print_arguments(args)

	if args.check_input:
		tg.print_args(args)
		print('Done checking input.')

	else:
		### set up TIGAR_Job class
		TIGAR_Job = Job.setup_job(args)

		### do job
		try:
			### do parallel jobs
			pool = multiprocessing.Pool(args.thread)
			pool.imap(TIGAR_Job.thread_process, [num for num in range(TIGAR_Job.n_targets)])
			pool.close()
			pool.join()
			print('Done.')

			### do any post-parallel jobs output file clean-up
			Job.clean_output(args, TIGAR_Job)

		finally:
			### other clean-up (currently only used by DPR training)
			Job.epilogue(args)

	### time calculation
	elapsed_sec = time() - start_time
	elapsed_time = tg.format_elapsed_time(elapsed_sec)
	print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

	### end logging
	sys.stdout.close()
