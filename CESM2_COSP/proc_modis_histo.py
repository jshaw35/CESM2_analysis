'''
Code by Jonah Shaw on 06/16/2021
Processing script to calculate cloud fraction histograms from 
MCD06COSP_M3_MODIS cloud products.

From Robert Pincus:
Thanks for the nudge. I enjoyed your talk yesterday at the AMWG meeting. I’ve cc’d Paul Hubanks (NASA) who put together the new data set.

As you’re likely aware there are two sets of pixel counts, one related to Cloud_Retrieval_Fraction_Total and the other to Cloud_Mask_Fraction. You should normalize the by the first. I include some notes based on an exchange between myself, Mark Zelinka (LLNL)  and Kerry Meyer (NASA) below. Let me know if there’s an remaining confusion.

"For this, I think your Option #1 is correct, using the sum of the non-PCL and PCL Cloud_Optical_Thickness joint histograms (vs CTP) as the numerator, and Cloud_Retrieval_Fraction_Total Pixel_Counts as the denominator. Because the cloud mask and the cloud optical properties retrievals define daytime slightly differently (using different thresholds on solar zenith angle), your denominator in Option #2 is a slightly larger “daytime” pixel population particularly near the terminator, and thus will yield fractions somewhat smaller than they should be.

Of course, using the Cloud_Optical_Thickness joint histograms will only yield the successful retrieval fractions that, because they do not include pixels for which the optical property retrievals fail, are inherently smaller than the cloud fraction derived from the cloud mask. But I assume from your past experience with our COSP products that you’re well aware of that.

Ideally, you could just use the Cloud_Retrieval_Fraction_Total mean itself without any manipulation, but that includes only our non-PCL [partly-cloudy and hence unreliable] pixels.”
''' 

from imports import *

def open_groups(filename):
    
    # Optical thickness / pressure histogram for cloudy retrievals
    taupres_histo = xr.open_dataset(filename,group='Cloud_Optical_Thickness_Total',engine='netcdf4')
    
    # Optical thickness / pressure histogram for partly cloudy retrievals
    taupres_histo_pcl = xr.open_dataset(filename,group='Cloud_Optical_Thickness_PCL_Total',engine='netcdf4')

    # Total counts for satellite retrievals (normalization)
    cld_counts = xr.open_dataset(filename,group='Cloud_Retrieval_Fraction_Total',engine='netcdf4')

    histo_normed = normalize_histo(taupres_histo, taupres_histo_pcl, cld_counts)
    
def normalize_histo(taupres_histo, taupres_histo_pcl, cld_counts):
    
    histo_counts_all = taupres_histo + taupres_histo_pcl # adds histograms for both cloudy and partly cloudy scenes
    
    histo_normed = histo_counts_all / cld_counts
    
    return histo_normed
    
#     cot_counts = cot_climo['JHisto_vs_Cloud_Top_Pressure'].sum(dim=['jhisto_cloud_optical_thickness_total_7','jhisto_cloud_top_pressure_7'])
#     cotpcl_counts = cotpcl_climo['JHisto_vs_Cloud_Top_Pressure'].sum(dim=['jhisto_cloud_optical_thickness_pcl_total_7','jhisto_cloud_top_pressure_7'])
    

    
    
# %%
if __name__ == '__main__':
    import sys

    
    
    args = sys.argv
    path_casefiles = Path(args[1])
    casename = path_casefiles.stem
    cs = CaseSetup(path_casefiles, casename)
    cs.create_case_all_tasks()
# %%