#!/usr/bin/env ruby
# 
#  preproc01.rb
#  
#  This script preprocessing rest and task-based fMRI data
#
#  Created by Zarrar Shehzad on 2012-09-24.
# 

# require 'pry'
# binding.pry

require 'pathname'
SCRIPTDIR   = Pathname.new(__FILE__).realpath.dirname.dirname
SCRIPTNAME  = Pathname.new(__FILE__).basename.sub_ext("")

# add lib directory to ruby path
$: << SCRIPTDIR + "lib" # will be scriptdir/lib

require 'config.rb'       # globalish variables
require 'for_commands.rb' # provides various function such as 'run'
require 'colorize'        # allows adding color to output
require 'erb'             # for interpreting erb to create report pages
require 'trollop'

# Process command-line inputs
p = Trollop::Parser.new do
  banner "Usage: #{File.basename($0)} -o orig-name -n scan -s sub1 ... subN\n"
  opt :orig, "Original dicom directory name piece", :type => :string, :required => true
  opt :name, "Name of scan (msit or rest)", :type => :string, :required => true
  opt :subjects, "Which subjects to process", :type => :strings, :required => true
end
opts = Trollop::with_standard_exception_handling p do
  raise Trollop::HelpNeeded if ARGV.empty? # show help screen
  p.parse ARGV
end

# Gather inputs
origsearch  = opts[:orig]
scan        = opts[:name]
subjects    = opts[:subjects]
ACQ         = "seq+z"

# html output    
layout_file       = SCRIPTDIR + "etc/layout.html.erb"
body_file         = SCRIPTDIR + "etc/01_preprocessing/#{SCRIPTNAME}.html.erb"
report_file       = "#{@qadir}/01_PreProcessed_#{SCRIPTNAME}_#{scan}.html"
@body             = ""

# Loop through each subject
subjects.each do |subject|
  puts "= Subject: #{subject} \n".white.on_blue
    
  puts "\n== Setting input variables".magenta
  in_subdir     = Dir["#{@origdir}/#{subject}*"].first
  in_rundirs    = Dir["#{in_subdir}/*#{origsearch}*"].sort
  runs          = 1..in_rundirs.length
  
  if in_rundirs.length == 0
    puts "ERROR: no runs found for search of #{in_subdir}/*#{origsearch}*".red
    next
  end
  
  puts "\n== Checking inputs".magenta
  next if any_inputs_dont_exist_including in_subdir
  
  puts "\n== Setting output variables".magenta
  out_subdir    = "#{@preprocdir}/#{subject}"
    
  puts "\n== Creating output directories (if needed)".magenta
  Dir.mkdir out_subdir if not File.directory? out_subdir
  
  runs.each_with_index do |run, i|
    puts "\n== Run #{run}".white.on_blue
    
    puts "\n=== Setting input variables".magenta
    in_rundir           = in_rundirs[i]
    
    puts "\n== Checking inputs".magenta
    next if any_inputs_dont_exist_including in_rundir    
    
    puts "\n=== Setting output variables".magenta
    out_scandir         = "#{out_subdir}/#{scan}"
    out_rundir          = "#{out_scandir}/run_%02d" % run
    orig                = "#{out_rundir}/orig.nii.gz"
    ppdir               = "#{out_rundir}/01_preprocess"
    mcref               = "#{out_subdir}/#{scan}/run_01/func_ref.nii.gz"  # same across runs
    motion              = "#{out_rundir}/motion.1D"
    motion_pic          = "#{out_rundir}/motion.png"
    abs_motion          = "#{out_rundir}/motion_absolute.1D"
    rel_motion          = "#{out_rundir}/motion_relative.1D"
    summary_motion      = "#{out_rundir}/motion_summary.txt"
    brain_mask          = "#{out_rundir}/func_mask.nii.gz"
    brain_axial_pic     = "#{out_rundir}/func_brain_mask_axial.png"
    brain_sagittal_pic  = "#{out_rundir}/func_brain_mask_sagittal.png"
    brain               = "#{out_rundir}/func_brain.nii.gz"
    mean                = "#{out_rundir}/func_mean.nii.gz"
    
    puts "\n== Saving contents for report page".magenta
    text      = File.open(body_file).read
    erbified  = ERB.new(text).result(binding)
    @body    += "\n #{erbified} \n"
    
    puts "\n== Checking outputs".magenta
    next if all_outputs_exist_including mcref, motion, abs_motion, rel_motion, 
                                        brain, brain_mask, mean, orig
    
    puts "\n== Creating output directories (if needed)".magenta
    Dir.mkdir out_scandir if not File.directory? out_scandir
    Dir.mkdir out_rundir if not File.directory? out_rundir
    Dir.mkdir ppdir if not File.directory? ppdir
    
    puts "\n== Converting from DICOM to NIFTI".magenta
    
    # NZ = number of images per volume assumes that files are 
    # named starting at 00001.dcm
    cmd = "dicom_hdr #{in_rundir}/00001.dcm | \
            grep '0019 100a' | \
            awk '{print $8}'"
    puts "#{cmd}".cyan
    NZ=`#{cmd}`.strip
    
    # TR = the repetition time (sampling frequency) of the
    # sequence
    cmd = "dicom_hdr #{in_rundir}/00001.dcm | \
            grep '0018 0080' | \
            awk '{print $9}' | sed 's/Time[/]*//g'"
    puts "#{cmd}".cyan
    TRms=`#{cmd}`.strip
    
    # sometimes afni prefers TR to be in seconds, so convert
    cmd = "echo 'scale=3;#{TRms}/1000' | bc -q 2>/dev/null"
    puts "#{cmd}".cyan
    TR=`#{cmd}`.strip
    
    # NT = the number of files, hopefully the number of TRs?
    cmd = "ls -1 #{in_rundir}/*.dcm | wc -l | tr '\n' ' '"
    puts "#{cmd}".cyan
    NT=`#{cmd}`.strip
    
    puts "?? epi(NT,NZ,TR) = #{NT},#{NZ},#{TR} ".light_magenta
    
    run "to3d -epan -session #{File.dirname(orig)} \
          -prefix #{File.basename(orig)} \
          -time:zt #{NZ} #{NT} #{TRms} '#{ACQ}' -assume_dicom_mosaic \
          #{in_rundir}/*.dcm"
    
    #puts "\n=== Excluding first 4 time points".magenta
    #run "3dcalc -a #{original}'[4..$]' -expr 'a' \
    #      -prefix #{ppdir}/01_exclude_tpts.nii.gz"
    
    puts "\n=== Performing slice time correction".magenta
    run "3dTshift -TR #{@TR} -slice #{@nslices} -tpattern #{@slice_pattern} \
          -prefix #{ppdir}/01_slice_time.nii.gz #{orig}"
    
    puts "\n=== Deobliquing to be AFNI friendly".magenta
    run "3dcopy #{ppdir}/01_slice_time.nii.gz #{ppdir}/02_deoblique.nii.gz"
    run "3drefit -deoblique #{ppdir}/02_deoblique.nii.gz"
    
    puts "\n=== Reorienting to be FSL friendly".magenta
    run "3dresample -inset #{ppdir}/02_deoblique.nii.gz -orient RPI \
          -prefix #{ppdir}/03_reorient.nii.gz"
    
    if i == 0
      run "3dcalc -a #{ppdir}/03_reorient.nii.gz'[0]' -expr 'a' -prefix #{mcref}"
    else
      run "ln -s #{mcref} #{out_rundir}/"
    end
    
    puts "\n=== Motion correcting".magenta
    puts "=== using the 1st image of run #1 as the reference".magenta
    run "3dvolreg -Fourier -prefix #{ppdir}/04_motion_correct.nii.gz \
          -base #{mcref} -1Dfile #{motion} #{ppdir}/03_reorient.nii.gz"
    
    puts "\n=== Calculating absolute and relative motion".magenta
    run "1d_tool.py -infile #{motion} -collapse_cols euclidean_norm \
          -overwrite -write #{abs_motion}"
    run "1d_tool.py -infile #{motion} -collapse_cols euclidean_norm \
          -derivative -overwrite -write #{rel_motion}"
    
    puts "\n=== Saving mean and max of relative motion".magenta
    rel_motion_contents = File.read(rel_motion).split("\n")
    relmean = rel_motion_contents.reduce(:+).to_f / rel_motion_contents.size
    relmax = rel_motion_contents.max
    File.open(summary_motion, 'w') {|file| file.write("# mean\tmax\n#{relmean}\t#{relmax}")}
    
    puts "\n=== Making pretty picture of motion".magenta
    run "fsl_tsplot -i #{abs_motion},#{rel_motion} -o #{motion_pic} \
          -t 'Estimated displacement (mm)' -a abs,rel"
    
    puts "\n=== Generating brain mask".magenta
    run "3dAutomask -dilate 1 -prefix #{brain_mask} #{ppdir}/04_motion_correct.nii.gz"
    
    puts "\n=== Creating pretty pictures".magenta
    run "slicer.py -w 5 -l 4 -s axial --overlay #{brain_mask} 1 1 -t #{mcref} #{brain_axial_pic}"
    run "slicer.py -w 5 -l 4 -s sagittal --overlay #{brain_mask} 1 1 -t #{mcref} #{brain_sagittal_pic}"
    
    puts "\n=== Applying mask to get only the brain".magenta
    run "3dcalc -a #{ppdir}/04_motion_correct.nii.gz -b #{brain_mask} \
          -expr 'a*ispositive(b)' -prefix #{brain}"
    
    puts "\n=== Creating average EPI".magenta
    run "3dTstat -mean -prefix #{mean} #{brain}"
  end
end

@title          = "Functional Preprocessing"
@nav_title      = @title
@dropdown_title = "Subjects"
@dropdown_elems = subjects
@foundation     = SCRIPTDIR + "lib/foundation"

puts "\n= Compiling and writing report page to %s".magenta % report_file
text      = File.open(layout_file).read
erbified  = ERB.new(text).result(binding)
File.open(report_file, 'w') { |file| file.write(erbified) }

