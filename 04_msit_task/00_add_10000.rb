#!/usr/bin/env ruby

# This script simply adds 10,000 to any data
# needed for inputing to FSL FEAT

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
  banner "Usage: #{File.basename($0)} -p ...\n"
  opt :paths, "File path to process", :type => :strings, :required => true
end
opts = Trollop::with_standard_exception_handling p do
  raise Trollop::HelpNeeded if ARGV.empty? # show help screen
  p.parse ARGV
end

paths = opts[:paths]

paths.each do |input|
  input = File.realpath(input)
  mask = "#{input.rmext}_mask.nii.gz"
  output = "#{input.rmext}_forfsl.nii.gz"
  
  run "fslmaths #{input} -abs -Tmin -bin #{mask}"
  run "fslmaths #{input} -add 100 -mas #{mask} #{output}"
  
  puts ""
end

