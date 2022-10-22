#!/bin/bsh

# Put something in this script that you want to happen when a new report has been rendered
# Remove the underscore prefix (_) from this file, and the file will be executed.
# You are given two arguments:
#  - 1: The title of the batch
#  - 2: The (absolute?) path to the report html file.


# This example sends an email
#if [ -z ${ASSCOM2_EMAIL+x} ]; then 
#    mail -s "ac2 report $1" -a $2 youreamailaddress@example.com <<< "The assemblycomparator2 report for project $1 has been rendered and is attached in this email."
#fi
