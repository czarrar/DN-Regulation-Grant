#!/usr/bin/env ruby

require 'pathname'
SCRIPTDIR   = Pathname.new(__FILE__).realpath.dirname.dirname
SCRIPTNAME  = Pathname.new(__FILE__).basename.sub_ext("")

# add lib directory to ruby path
$: << SCRIPTDIR + "lib" # will be scriptdir/lib


require 'config.rb'       # globalish variables
require 'for_commands.rb' # provides various function such as 'run'
require 'colorize'        # allows adding color to output
require 'trollop'         # command-line option parser

# Process command-line inputs
p = Trollop::Parser.new do
  banner "Usage: #{File.basename($0)} (-e 2) -n scan -r 1 .. N -s sub1 ... subN\n"
  opt :name, "Name of scan (movie or rest)", :type => :string, :required => true
  opt :runs, "Which runs to process", :type => :ints, :required => true
  opt :subjects, "Which subjects to process", :type => :strings, :required => true
end
opts = Trollop::with_standard_exception_handling p do
  raise Trollop::HelpNeeded if ARGV.empty? # show help screen
  p.parse ARGV
end

# Gather inputs
scan        = opts[:name]
runs        = opts[:runs]
subjects    = opts[:subjects]

roidir = "#{@projdir}/rois"
templates = "#{roidir}/PNAS_Smith09_rsn10_resample.nii.gz"
mask = "#{roidir}/#{scan}_mask.nii.gz"

subjects.each do |subject|
  
  puts "\n= Subject #{subject}".white.on_blue
    
  puts "\n== Setting input variables".magenta
  in_subdir     = "#{@preprocdir}/#{subject}"
  in_funcdir    = "#{in_subdir}/#{scan}"

  puts "\n== Setting output variables".magenta
  out_subdir    = "#{@analsubjdir}/#{subject}"
  out_funcdir   = "#{out_subdir}/#{scan}"
  
  Dir.mkdir out_subdir if not File.directory? out_subdir
  Dir.mkdir out_funcdir if not File.directory? out_funcdir
  
  runs.each do |run|
    
    puts "\n== Run #{run}".white.on_blue
    
    puts "\n=== Setting input variables".magenta
    in_rundir   = "#{in_funcdir}/run_%02d" % run
    func        = "#{in_rundir}/func_denoise+smooth2standard.nii.gz"
    
    puts "\nChecking inputs".magenta
    next if any_inputs_dont_exist_including func
    
    puts "\n=== Setting output variables".magenta
    out_rundir  = "#{out_funcdir}/run_%02d" % run
    tcs         = "#{out_rundir}/rsn10.1D"
    fcmaps      = "#{out_rundir}/rsn10.nii.gz"
    
    puts "\n=== Checking outputs".magenta
    next if all_outputs_exist_including tcs, fcmaps
    
    Dir.mkdir out_rundir if not File.directory? out_rundir
    
    puts "@@@ Extracting RSN TCs".magenta
    run "fsl_glm -i #{func} \
          -d #{templates} \
          -m #{mask} \
          -o #{tcs} \
          --demean"
    
    puts "@@@ Calculating RSN FS Maps".magenta
    run "fsl_glm -i #{func} \
          -d #{tcs} \
          -m #{mask} \
          -o #{fcmaps} \
          --demean"
  end
end
