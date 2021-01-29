#!/bin/sh
######################################################
#
# EXAMPLE SERIAL SCRIPT FOR SGE
# Modified by Basement Supercomputing 1/2/2006 DJE
# To use, change "MPICH_JOB", "NUMBER_OF_CPUS"
# and "MPICH_PROGRAM_NAME" to real values.
#
# Modified for NJIT - AM 25Jan07
#
######################################################
# Use private queue
#$ -q kenahn
#######################
# Job Name
#######################
#$ -N Dynamic_Initial_260_1
#
#########################################
# Send mail when jobs starts and finishes
#########################################
#$ -M lz242@njit.edu 
#$ -m be
###############################
# Use current working directory
###############################
#$ -cwd
#
###############################
# Join stdout and stderr
###############################
#$ -j y
#
###############################
# Declare job re-runnable
###############################
#$ -r y
###########################################################
# Run job through bash shell
###########################################################
#$ -S /bin/bash
##############################################################
# Full command path of your program
# --- Use full pathname to your executable ---
##############################################################
./main*