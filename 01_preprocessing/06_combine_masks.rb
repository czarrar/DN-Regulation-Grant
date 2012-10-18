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
end
opts = Trollop::with_standard_exception_handling p do
  raise Trollop::HelpNeeded if ARGV.empty? # show help screen
  p.parse ARGV
end

# Gather inputs
scan        = opts[:name]

puts "Merging masks together".magenta
infnames = Dir.glob("#{@preprocdir}/CCD*/#{scan}/run_*/func_mask2standard.nii.gz")
outfname = "#{@projdir}/rois/#{scan}_masks_all_subjects.nii.gz"
run "fslmerge -t #{outfname} #{infnames.join(' ')}"

puts "Getting mask with voxels present in all subjects".magenta
infname = outfname
outfname = "#{@projdir}/rois/#{scan}_mask.nii.gz"
run "fslmaths #{infname} -Tmin #{outfname}"
