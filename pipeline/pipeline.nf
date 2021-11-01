#!/usr/bin/env ms2rescore-immunopeptidomics-manuscript

/*
 * Define the (standard) parameters
 */

params.raw_file_location = false
params.mgf_store = false
params.raw_store = false
params.identifier = "no_identifier"
params.id_file = false
params.search_engine = false
params.config = "NO_FILE"
params.output_location = launchDir


output_location = file("$params.output_location/$params.identifier")
output_location.mkdir()
opt_config = file(params.config)

log.info """\
 ====================================================
 I M M U N O P E P T I D O M I C S    P I P E L I N E
 ====================================================
 Database identifier (PRIDE/massIVE)    : ${params.identifier}
 Raw Files Location                     : ${params.raw_file_location}
 Raw Files store Location               : ${params.raw_store}
 mgf Files store Location               : ${params.mgf_store}
 Peptide identification file            : ${params.id_file}
 Search engine used                     : ${params.search_engine}
 Output location                        : ${params.output_location}
 """


if (params.raw_file_location != false ){
    rawfiles = Channel.fromPath(params.raw_file_location)
}

else if (params.identifier =~ /^[P][X][D][0-9]{6}$/){

process Download_PRIDE_files{
    errorStrategy 'retry'
    if (params.raw_store == true){
        publishDir "$output_location/raw", mode: 'copy'
    }

    output:
    file '*.{raw,RAW,Raw}' into rawfiles

    script:
    """
    python ../immuno_ms2rescore_tools/download_pride_project.py $params.identifier -f raw
    """
}
}

else if (params.identifier =~ /^[M][S][V][0-9]{9}$/){

process Donwload_massIVE_files{
    errorStrategy 'retry'
    if (params.raw_store == true){
        publishDir "$output_location/raw", mode: 'copy'
    }

    output:
    file '*.{raw,RAW,Raw}' into rawfiles

    script:
    """
    python ../immuno_ms2rescore_tools/download_massive_project.py $params.identifier -f raw
    """
}
}

else {
    println "\n Invalid parameter(s): either give valid pxd/MSV identifier or raw file location"
    return;
}

process ThermoRawFileParser{

    stageInMode "copy"
    if (params.mgf_store == true){
       publishDir output_location/"mgf", mode: 'copy'
    }


    input:
    path raw_file from rawfiles.flatten()

    output:
    file "*.mgf" into mgffiles

    script:
    """
    thermorawfileparser --input=$raw_file -o=./ --format=0
    """
}

process IdFileParser{

    input:
    file idfile from Channel.fromPath(params.id_file, type: "file")
    file config from opt_config

    output:
    file "*.peprec" into peprec

    script:
    def filter = config.name != 'NO_FILE' ? "--config $config" : ""
    """
    python ../immuno_ms2rescore_tools/id_file_parser.py --id_file $idfile --search_engine $params.search_engine --output_filename $params.identifier
    """
}

process ConcatPeprecFiles{

    input:
    path peprec from peprec.collectFile()

    output:
    file peprec into final_peprec

    script:
    """
    sed -i '2,\${/^spec_id/d;}' $peprec
    """
}

process CreateSpectralLibary{
    publishDir output_location, mode: 'move'
    echo true

    input:
    path peprec from final_peprec
    path mgf_file from mgffiles.collect()

    output:
    file "*.peprec"
    file "*.mgf"

    script:
    """
    python ../immuno_ms2rescore_tools/spectral_library.py --peprec $peprec --mgf_folder . --identifier $params.identifier
    """
}

