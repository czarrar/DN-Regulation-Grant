#!/usr/bin/env ruby

# This script runs the first level stats to measure regions with
# signifiant modulation of BOLD signal as a result of the MSIT task
#
# It uses level1_template.fsf with feat
# and soft-links a reg directory at the end

require 'pathname'
SCRIPTDIR   = Pathname.new(__FILE__).realpath.dirname.dirname
SCRIPTNAME  = Pathname.new(__FILE__).basename.sub_ext("")

# add lib directory to ruby path
$: << SCRIPTDIR + "lib" # will be scriptdir/lib

require 'config.rb'       # globalish variables
require 'for_commands.rb' # provides various function such as 'run'
require 'colorize'        # allows adding color to output
require 'trollop'

# Process command-line inputs
p = Trollop::Parser.new do
  banner "Usage: #{File.basename($0)} -s sub1 ... subN\n"
  opt :subjects, "Which subjects to process", :type => :strings, :required => true
end
opts = Trollop::with_standard_exception_handling p do
  raise Trollop::HelpNeeded if ARGV.empty? # show help screen
  p.parse ARGV
end

subjects    = opts[:subjects]
runs = ["run_01", "run_02"]

out_basename = "task_hrf"
template = "level1_ccd_template.fsf"
subject_template = "CCD003"
run_template = "run_01"

# Loop through each subject
subjects.each do |subject|
  puts "= Subject: #{subject}".white.on_blue
  
  puts "= Setting IO variables".magenta
  in_subdir = "#{@preprocdir}/#{subject}"
  in_regdir = "#{in_subdir}/msit/reg"
  out_subdir = "#{@analsubjdir}/#{subject}"
  
  runs.each_with_index do |run, i|
    puts "== Run: #{i+1}".white.on_blue
    
    puts "== Setting IO variables".magenta
    out_rundir = "#{out_subdir}/msit/#{run}"
    out_fsf = "#{out_rundir}/#{out_basename}.fsf"
    out_featdir = "#{out_rundir}/#{out_basename}.feat"
    out_regdir = "#{out_featdir}/reg"
    
    puts "== Checking outputs".magenta
    next if all_outputs_exist_including out_fsf, out_featdir
    
    puts "Creating new fsf".magenta
    run "sed -e 's/#{subject_template}/#{subject}/g' -e 's/#{run_template}/#{run}/g' #{template} > #{out_fsf}"
    
    puts "Checking model setup".magenta
    run "feat_model #{out_rundir}/#{out_basename}"
    
    puts "Running model".magenta
    run "feat #{out_fsf}"
    
    puts "Soft-linking registration directory".magenta
    run "ln -s #{in_regdir} #{out_featdir}/"
    
    puts "Creating fake reg directory".magenta
    run "rm -r #{out_regdir}; echo 'good'"
    run "mkdir #{out_regdir}"
    run "cp /home2/data/PublicProgram/fsl-4.1.9//etc/flirtsch/ident.mat #{out_regdir}/example_func2standard.mat"
    run "cp /home2/data/PublicProgram/fsl-4.1.9/etc/flirtsch/ident.mat #{out_regdir}/highres2standard.mat"
    run "ln -sf /home2/data/PublicProgram/fsl-4.1.9//data/standard/MNI152_T1_4mm_brain.nii.gz #{out_regdir}/highres.nii.gz"
    run "ln -sf /home2/data/PublicProgram/fsl-4.1.9//data/standard/MNI152_T1_4mm_brain.nii.gz #{out_regdir}/standard.nii.gz"
    
  end
   
end