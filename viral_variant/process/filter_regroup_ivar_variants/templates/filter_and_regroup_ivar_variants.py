#!/usr/bin/env python3


import pandas as pd
import os



## Amend/Correct TOTAL_DP and ALT_FREQ from ivar for position with indels
def correct_variant_tsv(tsv_file):
    df = pd.read_csv(tsv_file, sep='\\t')

    grouped = df.groupby(['REGION', 'POS'])
    
    for (region, pos), group in grouped:        
        # correct only for coord with at least one indel
        with_indel = any(len(alt) > 1 for alt in group['ALT'])
        if with_indel:
            # recompute 'REF_DP' considering all lines (variant) on this coord
            total_dp = group['TOTAL_DP'].iloc[0]  # 'TOTAL_DP' is equal for all variant at this coord
            total_alt_dp = group['ALT_DP'].sum()
            corrected_ref_dp = total_dp - total_alt_dp
            df.loc[group.index, 'REF_DP'] = corrected_ref_dp

            # recompute 'REF_RV'
            total_dp_rv = group['REF_RV'].iloc[0]  # (TOTAL_DP == REF_DP) ==> (TOTAL_RV == REF_RV), can be false ??
            total_alt_dp_rv = group['ALT_RV'].sum()
            corrected_ref_dp_rv = total_dp_rv - total_alt_dp_rv
            df.loc[group.index, 'REF_RV'] = corrected_ref_dp_rv
    
    # save file    
    output_file = os.path.basename(tsv_file.replace('_raw', '').replace('.tsv', '')) + '_corrected.tsv'
    df.to_csv(output_file, sep='\\t', index=False)

    return df



def filter_variant_tsv(df, min_dp, ref_dp_ratio_max):
    filtered_df = df[(df['TOTAL_DP'] > min_dp) & 
                     (df['REF_DP'] / df['TOTAL_DP'] < ref_dp_ratio_max)]    
    return filtered_df



def summarize_pos(pos_list, out_prefix):
    # keep pos if filtered_in in at least one sample
    all_pos = pd.concat(pos_list).drop_duplicates()
    output_file = out_prefix + '_variable_pos.csv'
    all_pos.to_csv(output_file, sep='\\t', index=False)



def read_short_mpileup(short_mpileup_files):
    short_mpileup_dict = {}
    for short_mpileup_file in short_mpileup_files:
        sample_name = os.path.basename(short_mpileup_file.replace('_cols_1_4.mpileup', ''))
        mpileup_df = pd.read_csv(short_mpileup_file, sep='\\t', header=None, names=['REGION', 'POS', 'REF', 'TOTAL_DP'])
        short_mpileup_dict[sample_name] = mpileup_df
    return short_mpileup_dict



def filter_variants_by_batch(tsv_files, csv_positions, short_mpileup_dict, out_prefix):
    # import interest positions
    positions_df = pd.read_csv(csv_positions, sep='\\t')
    positions_set = set(zip(positions_df['REGION'], positions_df['POS']))
    
    output_dir = out_prefix + "_batchFiltered"
    os.makedirs(output_dir)
    common_columns = ['REGION', 'POS', 'REF', 'ALT', 'GFF_FEATURE', 'REF_CODON', 'REF_AA', 'ALT_CODON', 'ALT_AA', 'POS_AA']
    all_samples_filtered = []
    for tsv_file in tsv_files:
        sample_name = os.path.basename(tsv_file.replace('_raw', '').replace('.tsv', ''))
        ## import and filter
        tsv_df = pd.read_csv(tsv_file, sep='\\t', low_memory=False)     # 'low_memory=False' to avoid the following warning (which seems unfounded): '<stdin>:1: DtypeWarning: Columns (3,13,14,15,16,17,18) have mixed types. Specify dtype option on import or set low_memory=False.'. But it may have a significant impact on memory in the case of large genomes.
        print ("1111111111", tsv_df)
        filtered_df = tsv_df[tsv_df.apply(lambda row: (row['REGION'], row['POS']) in positions_set, axis=1)]
        print ("2222222222", filtered_df)

        # complete missing 'TOTAP_DP' about short_mpileup_dict
        if sample_name in short_mpileup_dict:
            mpileup_df = short_mpileup_dict[sample_name]
            filtered_df = filtered_df.merge(mpileup_df, on=['REGION', 'POS', 'REF'], how='left', suffixes=('', '_mpileup'))
            filtered_df['TOTAL_DP'] = filtered_df['TOTAL_DP'].fillna(filtered_df['TOTAL_DP_mpileup'])
            filtered_df.drop(columns=['TOTAL_DP_mpileup'], inplace=True)
        else:
            print(f"WARNING: no short_mpileup_cov file: '{sample_name}' ('TAOTAL_DP' not added for stable position).")
        print ("33333333333", filtered_df)

        # save in specific file complete filtered sample tsv
        filtered_output_path = os.path.join(output_dir, f"{sample_name}_batchFiltered.tsv")
        filtered_df.to_csv(filtered_output_path, sep='\\t', index=False)
        
        # reordone and rename specific columns using sample_name and append to list of df to merge
        variable_columns = [col for col in tsv_df.columns if col not in common_columns]
        filtered_df = filtered_df[common_columns + variable_columns]
        renamed_columns = {col: f"{col}({sample_name})" for col in variable_columns}
        filtered_df = filtered_df.rename(columns=renamed_columns)
        all_samples_filtered.append(filtered_df)

    # export filtered df of all sample on one unique file
    global_output_file = out_prefix + "_summary_all_iSNVs.tsv"
    global_df = all_samples_filtered[0]
    for df in all_samples_filtered[1:]:
        global_df = pd.merge(global_df, df, on=common_columns, how='outer')
    global_df.to_csv(global_output_file, sep='\\t', index=False)

    # TODO: for each pos in global_df, for each smpl, increade value of REF_DP, ALT_DP, ALT_FREQ and TOAL_DP using dictionay constructed about mpileup result.



def correct_and_filter_ivar_variant():
    tsv_files = "${ivar_tsv_files}".split()
    out_prefix = "${out_prefix}"
    min_dp = ${min_dp}
    ref_dp_ratio_max = ${ref_dp_ratio_max}

    all_filtered_pos = []
    for tsv_file in tsv_files:
        corrected_df = correct_variant_tsv(tsv_file)
        filtered_df = filter_variant_tsv(corrected_df, min_dp, ref_dp_ratio_max)
        output_file = os.path.basename(tsv_file.replace('_raw', '').replace('.tsv', '')) + '_filtered.tsv'
        filtered_df.to_csv(output_file, sep='\\t', index=False)
        all_filtered_pos.append(filtered_df[['REGION', 'POS']])

    summarize_pos(all_filtered_pos, out_prefix)



def regroup_ivar_variants():
    ### BE CAREFUL: correspondance between 'tsv_file' and 'short_mpileup_file' is based on basename, not robust, prefer use list index (e.g: list ordered following sample[meta.rank_in_batch], or simply follong meta.id (but not ideal order in final file))

    tsv_files_raw="${ivar_tsv_files}".split()
    tsv_files = [os.path.basename(tsv_file_raw.replace('_raw', '').replace('.tsv', '')) + '_corrected.tsv' for tsv_file_raw in tsv_files_raw]
    out_prefix = "${out_prefix}"
    csv_positions = "${out_prefix}_variable_pos.csv"
    short_mpileup_files = "${short_mpileup_cov_files}".split()
    
    short_mpileup_dict = read_short_mpileup(short_mpileup_files)
    filter_variants_by_batch(tsv_files, csv_positions, short_mpileup_dict, out_prefix)



if __name__ == "__main__":
    correct_and_filter_ivar_variant()
    regroup_ivar_variants()

