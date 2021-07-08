#!/usr/bin/perl 
# NOTE: The above path is machine dependent!!!

#------------------------------------------------------------------------
#  File: run_new_hfbrad_dme.pl
#
#  Authors:  R.J. Furnstahl  
#               furnstahl.1@osu.edu
#
#  Revision history: 
#   27-Apr-2008 --- original version, from run_hfbrad_dme.pl.  Allows a
#                    series of nuclei.
#   09-Jun-2008 --- added Coulomb switch

# This script runs the executable hfbrad_dme.x,
#  creating the appropriate input files and renaming the standard output files.
#  It is meant to run only a small number of "standard" nuclei as opposed
#  to large-scale runs. 

# Notes:
#  * Details on the standard HFBRAD inputs is given in 
#  * Everything on a line after a "#" is a comment.
#  * Make it executable with (do this if you get a "permission denied" error):
#       chmod +x run_new_hfbrad_dme.pl

# To do:
#  * make sure that dme_coefficients.inp gets linked (or checked)
#  * make it possible to put in multiple nuclei
#------------------------------------------------------------------------
#
#------------------------------------------------------------------------
# Prolog: Brief instructions in using this script
#  1. In the "setup" section, specify the nuclei to calculate in 
#      "nuclei_array" as a comma separated list.
#  2. Invoke the script as:
#      ./run_new_hfbrad_dme.pl <dme coefficient files>

use warnings;   # comment this out when not debugging the script

# version number
$version = "0.50  [09-Jun-2008]";

use File::Basename;   # functions to strip off the directory path of files
use Cwd;              # functions to find current directory 

#------------------------------------------------------------------------
# BEGIN: setup 
#------------------------------------------------------------------------

# From here until "end setup" you specify parameters of the calculation
@nuclei_array = ("O16","Ca40","Ca48","Ni56");  # comma separated list in quotes
#@nuclei_array = ("Ni56");  # comma separated list in quotes
#@nuclei_array = ("He4");  # comma separated list in quotes


$force = "Sn3l";      #  "SIII", "SLY4", "SKM*"
$force_string = $force; if ($force eq "SKM*") {$force_string = "SKMstar";};   
$mesh_points = 200;   # no. of r points           
$integ_step = 0.1;    # spacing of r points       
$it_max = 250;        # maximum no. of iterations 
$eps_energy = "1.e-9";  # energy tolerance
$max_delta = "1.e-7";   # maximum change tolerance
$boundary_condition = 0;
$xmu = 0.95;  # usually 0.65 or 0.7, but higher may be needed (0.8 works mostly)

$bogolyubov_p = "F";  # T if HFB, F if HF
$bogolyubov_n = "F";
$pairing_force = 3;
$regularization = "F"; # T if on, F if off

$coulomb_flag = "T";  # T for including Coulomb, F for no Coulomb

$dme_flag = "F";  # T for dme, F for regular Skyrme
$dme_type = "DME1";  # some label for diffent approximations
$kf_type = "LD";  #  LD = local density, TF = Thomas-Fermi, CB = Campi-Bouyssy

# $neutron and $proton now derived from hashes and @nuclei_array
#$neutron = 128;
#$proton = 82;
$j_max_n = 13;
$j_max_p = 13;


$densities = "T";
$quasiparticles = "F";


#------------------------------------------------------------------------
# END: setup
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# BEGIN: Definitions
#------------------------------------------------------------------------
%proton_hash = ( "He4" => 2,
                 "C12" => 6, 
		 "O16" => 8,
                 "Ca40" => 20,
		 "Ca48" => 20,
                 "Ni56" => 28,
                 "Pb208" => 82
                 );
	       
%neutron_hash = ( "He4" => 2,
                  "C12" => 6, 
		  "O16" => 8,
                  "Ca40" => 20,
		  "Ca48" => 28,
                  "Ni56" => 28,
                  "Pb208" => 128
                 );
#------------------------------------------------------------------------
# END: Definitions
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# BEGIN: Process command-line options
#------------------------------------------------------------------------

# set default options
$opt_version = 0;      # don't print version
$opt_help = 0;         # don't print help

# get the command-line options
use Getopt::Long;   # allow the extended processing of command line options
GetOptions("help", "version", "nodes=i", "ppn=i", "processors=i",
           "iterations=i", "compiler=s", "IO!");

# print out version if requested, then quit
if ($opt_version) {
  print "\nVersion number: $version\n";
  exit 1;
}

# print out help if requested, then quit
if ($opt_help) {
  print_help();
  exit 1;
}

# check for file names; quit if none specified
if (!$ARGV[0]) { 
  print "You must specify at least one dme input file.\n"; 
  print "Usage: run_new_hfbrad_dme.pl [--help] [--version] ";
  print " file(s)\n";
  exit 1;
}
  
#------------------------------------------------------------------------
# END: Process command-line options
#------------------------------------------------------------------------



#------------------------------------------------------------------------
# BEGIN: Generate the input file and run the code 
#------------------------------------------------------------------------

# Step through the dme input files
foreach $dme_input_file_full (@ARGV) {

  # put a link to this file
  print "Processing $dme_input_file_full ...\n";
  unlink("dme_coefficients.inp");
  symlink($dme_input_file_full,"dme_coefficients.inp");

  $dme_output_base = basename($dme_input_file_full);  # strip any directories
  $dme_output_base =~ s/\.out//;    # strip the ".out" ending

  foreach $nucleus (@nuclei_array) {  # step through the nuclei
    print "  Doing $nucleus ... \n";
    $proton = $proton_hash{$nucleus};
    $neutron = $neutron_hash{$nucleus};
    $file_np = sprintf("_N%1d_Z%1d",$neutron,$proton);

    $output_dir = "new_output/output_test_" . $dme_output_base . "_" . $force_string . $file_np;
    print "Output directory: $output_dir \n ";  # exit 1;

    if (-d $output_dir) {
      print "Output directory $output_dir already exists.  Files will be overwritten.\n";
    } else {
      mkdir($output_dir) || die "Can't make $output_dir directory.";
    }

    $file_log = "hfbrad" . "_" . $force_string . $file_np . ".log";
    $file_spe = "hfb" . "_" . $force_string . $file_np . ".spe";
    $file_neutron_density = "neutron" . "_" . $force_string . $file_np . ".dens";
    $file_proton_density = "proton" . "_" . $force_string . $file_np . ".dens";
    $file_mean_fields = "pot" . "_" . $force_string . $file_np . ".mf";
    $file_summary = "hfb" . "_" . $force_string . $file_np . ".summary"; 
    $file_input = "hfbrad" . "_" . $force_string . $file_np . ".inp"; 

    open(INFILE, ">hfb.input");

      print INFILE "&input \n";
      print INFILE "force = \"$force\", \n";
      print INFILE "mesh_points = $mesh_points, \n";
      print INFILE "integ_step = $integ_step, \n";
      print INFILE "it_max = $it_max, \n";
      print INFILE "eps_energy = $eps_energy, \n";
      print INFILE "max_delta = $max_delta, \n";
      print INFILE "boundary_condition = $boundary_condition, \n";
      print INFILE "xmu = $xmu, \n";
      print INFILE "bogolyubov = $bogolyubov_p, $bogolyubov_n, \n";
      print INFILE "pairing_force = $pairing_force, \n";
      print INFILE "regularization = $regularization \n";
      print INFILE "/ \n";

      print INFILE " \n";

      print INFILE "&dme_namelist \n";
      print INFILE "dme_flag = $dme_flag, \n";
      print INFILE "dme_type = $dme_type, \n";
      print INFILE "kf_type = $kf_type \n";
      print INFILE "/ \n";

      print INFILE " \n";

      print INFILE "&nucleus \n";
      print INFILE "neutron = $neutron, \n";
      print INFILE "proton = $proton, \n";
      print INFILE "j_max = $j_max_n, $j_max_p, \n";
      print INFILE "densities = $densities, \n";
      print INFILE "quasiparticles = $quasiparticles \n";
      print INFILE "/ \n";

    close(INFILE);

    # Run the code
    print "     Running HFBRAD with force $force";
    print " on nucleus N = $neutron, Z = $proton . . .\n";
    `./hfbrad_dme.x > $file_log`;

    # Rename output files
    print "\n *** All done! ***\n";
    print "Moving everything to $output_dir . . . \n";

    print "input file = $file_input \n";
    rename ("hfb.input", $output_dir . "/" . $file_input)
      || die "cannot rename input file";

    print "log file = $file_log \n";
    rename ($file_log, $output_dir . "/" . $file_log)
      || die "cannot rename log file";

    print "summary file = $file_summary \n";
    rename ("hfb.summary", $output_dir . "/" . $file_summary)
      || die "cannot rename summary file";

    print "spe file = $file_spe \n";
    rename ("hfb" . $file_np . ".spe", $output_dir . "/" . $file_spe)
      || die "cannot rename spe file";

    print "neutron density file = $file_neutron_density \n";
    rename ("neutron" . $file_np . ".dens", $output_dir . "/" . $file_neutron_density)
      || die "cannot rename neutron density file";

    print "proton density file = $file_proton_density \n";
    rename ("proton" . $file_np . ".dens", $output_dir . "/" . $file_proton_density)
      || die "cannot rename proton density file";

    print "Kohn-Sham potentials file = $file_mean_fields \n";
    rename ("pot" . $file_np . ".mf", $output_dir . "/" . $file_mean_fields)
      || die "cannot rename mean fields file";

    # Move wavefunction files to output directory  
    `/bin/mv wf* $output_dir`;  

  }   # close foreach $nucleus loop

}  # close foreach $dme_input_file_full loop 

#------------------------------------------------------------------------
# END: Run the code 
#------------------------------------------------------------------------
