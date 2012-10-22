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
require 'trollop'         # command-line parsing

# Process command-line inputs
p = Trollop::Parser.new do
  banner "Usage: #{File.basename($0)} -s sub1 ... subN\n"
  opt :subjects, "Which subjects to process", :type => :strings, :required => true
end
opts = Trollop::with_standard_exception_handling p do
  raise Trollop::HelpNeeded if ARGV.empty? # show help screen
  p.parse ARGV
end

# Gather inputs
subjects    = opts[:subjects]

# html output    
layout_file       = SCRIPTDIR + "etc/layout.html.erb"
body_file         = SCRIPTDIR + "etc/01_preprocessing/#{SCRIPTNAME}.html.erb"
report_file       = "#{@qadir}/01_PreProcessed_#{SCRIPTNAME}.html"
@body             = ""

# Loop through each subject
subjects.each do |subject|
  puts "\n= Subject: #{subject} \n".white.on_blue
  
  puts "\n== Setting input variables".magenta
  in_subdir       = Dir["#{@origdir}/#{subject}_*"].first
  in_anatdir      = "#{in_subdir}/DMN_FB_ANAT_2"
  
  puts "\n== Checking inputs".magenta
  next if any_inputs_dont_exist_including in_anatdir
  
  puts "\n== Setting output variables".magenta
  out_subdir      = "#{@preprocdir}/#{subject}"
  mridir          = "#{@freesurferdir}/#{subject}/mri"
  anatdir         = "#{out_subdir}/anat"
  orig            = "#{anatdir}/orig.nii.gz"
  head            = "#{anatdir}/head.nii.gz"
  brain           = "#{anatdir}/brain.nii.gz"
  brain_mask      = "#{anatdir}/brain_mask.nii.gz"  
  head_axial_pic          = "#{head.rmext}_axial_pic.png"
  head_sagittal_pic       = "#{head.rmext}_sagittal_pic.png"
  brain_mask_axial_pic    = "#{brain_mask.rmext}_axial_pic.png"
  brain_mask_sagittal_pic = "#{brain_mask.rmext}_sagittal_pic.png"
  
  puts "\n== Saving contents for report page".magenta
  text      = File.open(body_file).read
  erbified  = ERB.new(text).result(binding)
  @body    += "\n #{erbified} \n"
  
  puts "\n== Checking outputs".magenta
  next if all_outputs_exist_including head, brain, brain_mask
      
  puts "\n== Creating output directories (if needed)".magenta
  Dir.mkdir out_subdir if not File.directory? out_subdir
  Dir.mkdir anatdir if not File.directory? anatdir
    
  begin
    
    puts "\n== Converting DICOM to NIFTI".magenta
    run "to3d -anat -session #{File.dirname(orig)} -prefix #{File.basename(orig)} #{in_anatdir}/*.dcm"
    
    if File.exists? "#{mridir}/T1.mgz" and File.exists? "#{mridir}/brainmask.mgz"
      puts "\n== Freesurfer output already exists, skipping recon-all".red
    else      
      puts "\n== Using freesurfer for skull stripping".magenta
      run "recon-all -i #{orig} -autorecon1 -s #{subject} -sd #{@freesurferdir}"
    end
    
    puts "\n== Converting freesurfer output to nifti format".magenta
    run "mri_convert #{mridir}/orig.mgz #{anatdir}/tmp_head.nii.gz"
    run "mri_convert #{mridir}/brainmask.mgz #{anatdir}/tmp_brain.nii.gz"
    
    puts "\n== Reorienting head and to be FSL friendly".magenta
    run "3dresample -orient RPI -inset #{anatdir}/tmp_head.nii.gz -prefix #{head}"
    run "3dresample -orient RPI -inset #{anatdir}/tmp_brain.nii.gz -prefix #{anatdir}/tmp_brain_reorient.nii.gz"
    
    puts "\n== Generating brain mask and brain".magenta
    run "3dcalc -a #{anatdir}/tmp_brain_reorient.nii.gz -expr 'step(a)' -prefix #{brain_mask}"
    run "3dcalc -a #{head} -b #{brain_mask} -expr 'a*step(b)' -prefix #{brain}" 
    
    puts "\n== Creating pretty pictures".magenta
    run "slicer.py -w 5 -l 4 -s axial #{head} #{head_axial_pic}"
    run "slicer.py -w 5 -l 4 -s sagittal #{head} #{head_sagittal_pic}"
    run "slicer.py -w 5 -l 4 -s axial --overlay #{brain_mask} 1 1 -t #{head} #{brain_mask_axial_pic}"
    run "slicer.py -w 5 -l 4 -s sagittal --overlay #{brain_mask} 1 1 -t #{head} #{brain_mask_sagittal_pic}"
    
  ensure
    
    puts "\n== Removing intermediate files".magenta
    run "rm -f #{anatdir}/tmp_*.nii.gz"
    
  end
  
end


@title          = "Anatomical Preprocessing"
@nav_title      = @title
@dropdown_title = "Subjects"
@dropdown_elems = subjects
@foundation     = SCRIPTDIR + "lib/foundation"

puts "\n= Compiling and writing report page to %s".magenta % report_file
text      = File.open(layout_file).read
erbified  = ERB.new(text).result(binding)
File.open(report_file, 'w') { |file| file.write(erbified) }
