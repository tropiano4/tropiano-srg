#!/usr/bin/env python3

#------------------------------------------------------------------------------
# File: test.py
#
# Author:   A. J. Tropiano (tropiano.4@osu.edu)
# Date:     May 1, 2019
# 
# The purpose of this script is to test out codes. Generally, something is
# tested here, which evolves into several tests/checks. In the case of an
# extensive script, one should save the script as a seperate file with an
# extension _testv#.py where v# corresponds to the version number. For example,
# momentum_projection_operator_testv1.py. Use the revision history below to
# document when and why these files are created.
#
# Revision history:
#   08/28/19 --- Testing the deuteron momentum distribution using projection
#                operators from operators.py
#   08/30/19 --- Testing how to use *args in Python functions. This may be
#                useful for plotting codes. Created
#                momentum_projection_operator_testv1.py in Old_codes based off
#                last tests in this script.
#
#------------------------------------------------------------------------------


def func(channel, *args):

    print(channel)
    for arg in args:
        print('kvnn = %d'%arg[0])
        print('kmax = %.1f'%arg[1])
        print('kmid = %.1f'%arg[2])
        print('generator = %s'%arg[3])
        if arg[3] == 'Block-diag':
            print('Lambda_bd = %.2f'%arg[4])

channel = '3S1'

em_n3lo_wegner = (10, 30.0, 4.0, 'Wegner')
rke_n3lo_wegner = (106, 8.0, 2.0, 'Wegner')
gez_n2lo_wegner = (222, 10.0, 2.0, 'Block-diag', 2.00)

func(channel, em_n3lo_wegner, rke_n3lo_wegner, gez_n2lo_wegner)