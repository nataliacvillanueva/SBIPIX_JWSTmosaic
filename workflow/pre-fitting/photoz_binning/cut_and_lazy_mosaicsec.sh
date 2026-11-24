#!/bin/bash
read -p "Enter number of sections to split into: " tot_num_sections

for section in $(seq 0 $((tot_num_sections - 1))); do
    prefix_file_name="lazymosaic_sec${section}"

    echo "---------------------------------------"
    echo "Generating catalog for section ${section} of $((tot_num_sections - 1)) total sections."
    echo "Using filename prefix: ${prefix_file_name}"
    python makecatalog_general.py --output_fname_prefix "$prefix_file_name" --section_number "$section" --tot_sections "$tot_num_sections"

    gsed -i -E "/^[[:space:]]*\[io\][[:space:]]*$/,/^[[:space:]]*\[.*\][[:space:]]*$/ {
        /^[[:space:]]*output_file[[:space:]]*=/ {
        s|^([[:space:]]*output_file[[:space:]]*=[[:space:]]*)[^#[:space:]]+([[:space:]]*(#.*)?)$|\1\"./lazymosaic_catalogsections/${prefix_file_name}_outlazy.fits\"\2|
        }
    }" my_params.toml

    gsed -i -E "/^[[:space:]]*\[io\][[:space:]]*$/,/^[[:space:]]*\[.*\][[:space:]]*$/ {
    /^[[:space:]]*input_catalog[[:space:]]*=/ {
        s|^([[:space:]]*input_catalog[[:space:]]*=[[:space:]]*)[^#[:space:]]+([[:space:]]*(#.*)?)$|\1\"./lazymosaic_catalogsections/${prefix_file_name}_inlazy.fits\"\2|
    }
    }" my_params.toml

    echo "Running Lazy fit with updated parameter file..."

    lazy fit -p my_params.toml -t auto

    echo "Completed section ${section}."
done