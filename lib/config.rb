require 'for_commands.rb'

# Set paths
ENV['BASEDIR']  ||= "/home2/data"
@basedir          = ENV['BASEDIR']
@study            = "CCD"
@projdir          = "#{@basedir}/Projects/#{@study}"
@origdir          = "#{@projdir}/originals"
@qadir            = "#{@projdir}/qc"
@preprocdir       = "#{@projdir}/preprocessed"
@freesurferdir    = "#{@projdir}/freesurfer"
@roidir           = "#{@projdir}/rois"
@analdir          = "#{@projdir}/analysis"
@analsubjdir      = "#{@analdir}/subjects"
@priordir         = "/home2/data/PublicProgram/C-PAC/tissuepriors/2mm"

@TR               = 2.0
@nslices          = 17
@slice_pattern    = "seq+z"

exit 2 if any_inputs_dont_exist_including @origdir, @qadir, @preprocdir, 
                                          @freesurferdir
