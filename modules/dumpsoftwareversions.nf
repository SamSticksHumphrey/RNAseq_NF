//
//
//Custom module used to dump software versions within the nf-core pipeline template

process CUSTOM_DUMPSOFTWAREVERSIONS {
  publishDir = [
      path: { "${params.batch}" },
      mode: 'copy',
      pattern: '*_versions.yml'
  ]

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    script:
    template 'dumpsoftwareversions.py'



    stub:
    """
    touch software_versions.yml
    touch software_versions_mqc.yml
    touch versions.yml
    """

}
