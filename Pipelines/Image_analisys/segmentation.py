def segmentation(filename, params):

    import numpy as np
    from scipy import stats
    import time, os, sys
    from urllib.parse import urlparse
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import colorsys
    mpl.rcParams['figure.dpi'] = 300
    from cellpose import plot, models, utils, io
    from skimage import measure, filters
    from skimage.filters import threshold_otsu
    from skimage import restoration
    from skimage.segmentation import watershed
    from scipy import ndimage as ndi
    from skimage.feature import peak_local_max
    from skimage import img_as_ubyte, img_as_float
    from PIL import Image
    import csv
    from tqdm import tqdm
    from time import sleep
    
    # Nuclear signal
    nu_sig=1
    
    # Cytoplasmic signal
    cell_sig=0
    
    # Protein signal
    prot_sig=2
    
    # Nuclear image
    vmin_nu=0.1
    vmax_nu=99.7
    
    # Cell image
    vmin_cell=0.1
    vmax_cell=99
    
    # Protein images
    vmin_prot=0.1
    vmax_prot=99
    
    # Read file
    image=io.imread(filename)
    
    # Create file for analysis (channel1=cell_signal,channel2=nuclei,channel3=protein_signal)
    # Need to check if there it is a z-stack:
    # if there is no z --> (channel,x,y): dim=3
    # if there is z --> (channel,z,x,y): dim=4 --> perform max projection
    
        
    image_shp=np.array(image.shape).shape[0]
    if (image_shp==3):    
        temp=(image[cell_sig],image[nu_sig],image[prot_sig])
        img=np.asarray(temp)
    if (image_shp==4):
        image_proj=np.max(image,axis=0)
        temp=(image_proj[cell_sig],image_proj[nu_sig],image_proj[prot_sig])
        img=np.asarray(temp)
           
    # Extract channels
    cell_raw=img[0].copy()
    nuclei_raw=img[1].copy()
    protein_raw=img[2].copy()

    ################### CELLS SEGMENTATION ##############################
    
    # Gamma correction (I^gamma)
    cell_gamma=np.power(cell_raw,params.get('gamma_cell'))
    # Gaussian blur
    cell_blur=filters.gaussian(cell_gamma, sigma=params.get('gaussian_cell'))
    
    
    # Define CellPose Model
    model_cell = models.Cellpose(gpu=params.get('use_gpu'), model_type=params.get('model_type'))
    
    # Cellpose sementation using cell + nuclear channel or only cell channel
    if (params.get('nu_in_cell_seg') == True):
        
        # Channels definition
        ch = [1,2]
        
        # Gaussian blur with the same sigma as cells
        nuclei_blur_cell_seg=filters.gaussian(nuclei_raw, sigma=params.get('gaussian_cell'))
        
        # Image definition
        cell_blur=np.asarray((cell_blur,nuclei_blur_cell_seg))
    
        print("Segmentation of cytoplasm with CellPose using cell + nuclear signal...", flush=True)
    
        # Perform segmentation
        masks_cell, flows_cell, styles_cell, diams_cell = model_cell.eval(cell_blur, diameter=params.get('estimated_cell_diameter'), channels=ch, flow_threshold=params.get('flow_threshold'))
    
        print("Done!", flush=True)
    
    else:
        
        # Channels definition (no nuclei)
        ch = [0,0]
    
        print("Segmentation of cytoplasm with CellPose using cyto only...", flush=True)
    
        # Perform segmentation
        masks_cell, flows_cell, styles_cell, diams_cell = model_cell.eval(cell_blur, diameter=params.get('estimated_cell_diameter'), channels=ch, flow_threshold=params.get('flow_threshold'))
    
        print("Done!", flush=True)
    
    
    masks_cell_fil_2 = masks_cell.copy()
    print("Filtering cell area...", flush=True)
    props_cell_fil_2 = measure.regionprops(masks_cell_fil_2,cell_raw)
    masks_cell_fil_3=masks_cell_fil_2.copy()
    for i in range(0,len(props_cell_fil_2)):
        if (props_cell_fil_2[i].area<params.get('filter_area_cell')):
            masks_cell_fil_3[masks_cell_fil_3==props_cell_fil_2[i].label]=0
                
    # Final mask
    masks_cell_fil=masks_cell_fil_3.copy()
    
    io.save_masks(cell_gamma, 
              masks_cell_fil, 
              flows_cell, 
              filename, 
              channels=ch,
              png=True, # save masks as PNGs and save example image
              tif=False, # save masks as TIFFs
              save_txt=False, # save txt outlines for ImageJ
              save_flows=False, # save flows as TIFFs
              save_outlines=True, # save outlines as TIFFs 
              save_mpl=False # make matplotlib fig to view (WARNING: SLOW W/ LARGE IMAGES)
              )
    print("Done!", flush=True)
    ################### NUCLEI SEGMENTATION ##############################

    # Gamma correction (x^gamma)
    nuclei_gamma=np.power(nuclei_raw,params.get('gamma_nu'))
    # Gaussian blur
    nuclei_blur=filters.gaussian(nuclei_gamma, sigma=params.get('gaussian_nu'))
    
    print("Segmentation of Nuclei with Cellpose...", flush=True)
    # Define channels
    ch_nu=[0,0]
    # Define Cellpose nuclei model
    model_nu = models.Cellpose(gpu=params.get('use_gpu'), model_type=params.get('model_type_nucleus'))
    # Run evaluation for nuclei segmentation
    masks_nu, flows_nu, styles_nu, diams_nu = model_nu.eval(nuclei_blur, diameter=params.get('estimated_nucleus_diameter'), channels=ch_nu, flow_threshold=params.get('flow_threshold'))
    
    masks_nu_fil_2=masks_nu.copy()
    print("Done!", flush=True)
    
    print("Filtering nuclei area and eccentricity...", flush=True)
    props_nu_fil_2 = measure.regionprops(masks_nu_fil_2,nuclei_raw)
    masks_nu_fil_3=masks_nu_fil_2.copy()
    for i in range(0,len(props_nu_fil_2)):
        if (props_nu_fil_2[i].area<params.get('filter_area_nu')):
            masks_nu_fil_3[masks_nu_fil_3==props_nu_fil_2[i].label]=0
        else:
            if (props_nu_fil_2[i].eccentricity>params.get('filter_ecc')):
                masks_nu_fil_3[masks_nu_fil_3==props_nu_fil_2[i].label]=0
                
    # Final mask
    masks_nu_fil=masks_nu_fil_3.copy()
    
    print("Done!", flush=True)
    
    io.save_masks(nuclei_blur, 
                  masks_nu_fil, 
                  flows_nu, 
                  file_names = os.path.splitext(filename)[0]+'_nuclei',
                  channels=ch_nu,
                  png=True, # save masks as PNGs and save example image
                  tif=False, # save masks as TIFFs
                  save_txt=False, # save txt outlines for ImageJ
                  save_flows=False, # save flows as TIFFs
                  save_outlines=True, # save outlines as TIFFs 
                  save_mpl=False # make matplotlib fig to view (WARNING: SLOW W/ LARGE IMAGES)
                  )

    ################### FILTERING ############################## 

    # Create color map for plotting masks
    n=2**16
    h = (0,1)
    l = (.4,1)
    s =(.2,.8)
    h,l,s = np.random.uniform(*h,n), np.random.uniform(*l,n), np.random.uniform(*s,n)
    cols = np.stack([colorsys.hls_to_rgb(_h,_l,_s) for _h,_l,_s in zip(h,l,s)],axis=0)
    cols[0] = 0 
    mycmap = mpl.colors.ListedColormap(cols) 
    pixel_size = params.get('pixel_size')
    
    print("Filtering cells and nuclei...", flush=True) 
    # Filtering of nuclei
    props_nu_fil = measure.regionprops(masks_nu_fil,nuclei_raw)
    masks_nu_cell_fil_1 = masks_nu_fil.copy() 
    
    for i in tqdm(range(0,len(props_nu_fil))):
        #print(f'{i} of {len(props_nu_fil)}')
        
        masks_single_nu_in_cell=np.where(masks_nu_fil==props_nu_fil[i].label,masks_cell_fil,0)
        props_masks_single_nu_in_cell=measure.regionprops(masks_single_nu_in_cell,cell_raw)
        
        if (len(props_masks_single_nu_in_cell)>2):
            masks_nu_cell_fil_1[masks_nu_fil==props_nu_fil[i].label]=0
        
        
        if (len(props_masks_single_nu_in_cell)==2):
            masks_single_cell1=np.where(masks_cell_fil==props_masks_single_nu_in_cell[0].label,masks_cell_fil,0)
            props_masks_single_cell1=measure.regionprops(masks_single_cell1,cell_raw)
            masks_single_cell2=np.where(masks_cell_fil==props_masks_single_nu_in_cell[1].label,masks_cell_fil,0)
            props_masks_single_cell2=measure.regionprops(masks_single_cell2,cell_raw)
            if (props_masks_single_nu_in_cell[0].area>=props_masks_single_cell1[0].area or props_masks_single_nu_in_cell[1].area>=props_masks_single_cell2[0].area):
                 masks_nu_cell_fil_1[masks_nu_fil==props_nu_fil[i].label]=0
            else:
                if (min(props_masks_single_nu_in_cell[0].area,props_masks_single_nu_in_cell[1].area)/max(props_masks_single_nu_in_cell[0].area,props_masks_single_nu_in_cell[1].area)>0.8):
                    masks_nu_cell_fil_1[masks_nu_fil==props_nu_fil[i].label]=0
                else:
                    if (props_masks_single_nu_in_cell[0].area>props_masks_single_nu_in_cell[1].area):
                        masks_nu_cell_fil_1[masks_single_nu_in_cell==props_masks_single_nu_in_cell[1].label]=0
                    else:
                        masks_nu_cell_fil_1[masks_single_nu_in_cell==props_masks_single_nu_in_cell[0].label]=0
        
        
        if (len(props_masks_single_nu_in_cell)==1):
            if (props_masks_single_nu_in_cell[0].area<params.get('filter_area_nu')):
                masks_nu_cell_fil_1[masks_nu_fil==props_nu_fil[i].label]=0
            else:
                masks_single_cell=np.where(masks_cell_fil==props_masks_single_nu_in_cell[0].label,masks_cell_fil,0)
                props_masks_single_cell=measure.regionprops(masks_single_cell,cell_raw)
                if (props_masks_single_nu_in_cell[0].area>=props_masks_single_cell[0].area):
                    masks_nu_cell_fil_1[masks_nu_fil==props_nu_fil[i].label]=0
        
        
        if (len(props_masks_single_nu_in_cell)==0):
                masks_nu_cell_fil_1[masks_nu_fil==props_nu_fil[i].label]=0
                
    # Filtering of cells            
    props_cell = measure.regionprops(masks_cell_fil,cell_raw)

    for i in range(len(props_cell)):
        masks_single_cell=np.where(masks_cell_fil==props_cell[i].label,masks_cell_fil,0)
        masks_nu_in_single_cell=np.where(masks_single_cell!=0,masks_nu_cell_fil_1,0)
        props_nu_in_single_cell=measure.regionprops(masks_nu_in_single_cell,nuclei_raw)
        if (len(props_nu_in_single_cell)==0):
            masks_cell_fil[masks_cell_fil==props_cell[i].label]=0
           
    # Final masks generation
    #masks_cell_fil=masks_cell_fil_1.copy()
    masks_nu_cell_fil=masks_nu_cell_fil_1.copy()
    
    # Create masks of nuclei with the same index as cells
    masks_nu_in_cell=np.where(masks_nu_cell_fil!=0,masks_cell_fil,0)
    
    print("Done!", flush=True)
    
    # Create masks of cytoplasm where there is no nucleus
    masks_cyto_only=np.where(masks_nu_cell_fil==0,masks_cell_fil,0)
    
    # Calculate properties for protein in cell, nucleus, cytoplasm
    props_protein_in_cyto = measure.regionprops(masks_cyto_only,protein_raw)
    props_protein_in_cell = measure.regionprops(masks_cell_fil,protein_raw)
    props_protein_in_nu = measure.regionprops(masks_nu_in_cell,protein_raw)
    
    output_folder = params.get('output_directory')
    base_name = os.path.splitext(os.path.basename(filename))[0]
    output_dir = os.path.join(output_folder, f"{base_name}_output")
    
    finalmasks_dir = os.path.join(os.path.dirname(filename), f'{base_name}_filtered_masks')
    
    # Create output folder if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(finalmasks_dir):
        os.makedirs(finalmasks_dir)
    
    # 1. Nuclear Segmentation
    # Create Figure Object
    plt.figure(figsize=(6,6))
    # Plot Raw nuclei with selected percentiles 
    plt.imshow(nuclei_raw+cell_raw, vmin=np.percentile(nuclei_raw,min(vmin_nu,vmin_cell)), vmax=np.percentile(nuclei_raw,max(vmax_nu,vmax_cell)), cmap='gray')
    # Plot Nuclear segmentation masks
    plt.imshow(mycmap(masks_nu_in_cell), alpha=0.5)
    # Do not plot axis
    plt.axis("off") 
    # Set title
    plt.title("Nuclear Masks", fontsize=20, color='white')
    # Save figure
    plt.savefig(output_dir+'/nuclear_segmentation.tif',transparent=True)
    # Close figure
    plt.close()
    print("1. Exported nuclear masks on cell + nuclear signal!", flush=True)
    
    # 2. Cell Segmentation
    io.save_masks(cell_gamma, 
              masks_cell_fil, 
              flows_cell, 
              file_names = os.path.join(finalmasks_dir, os.path.basename(os.path.splitext(filename)[0])+'_final_cell_mask'),
              channels=ch,
              png=True, # save masks as PNGs and save example image
              tif=False, # save masks as TIFFs
              save_txt=False, # save txt outlines for ImageJ
              save_flows=False, # save flows as TIFFs
              save_outlines=True, # save outlines as TIFFs 
              save_mpl=False # make matplotlib fig to view (WARNING: SLOW W/ LARGE IMAGES)
              )
    # Create Figure Object
    plt.figure(figsize=(6,6))
    # Plot Raw cytoplasm with selected percentiles 
    plt.imshow(cell_raw+nuclei_raw, vmin=np.percentile(cell_raw,min(vmin_nu,vmin_cell)), vmax=np.percentile(cell_raw,max(vmax_nu,vmax_cell)), cmap='gray')
    # Plot cytoplasm segmentation masks
    plt.imshow(mycmap(masks_cell_fil), alpha=0.5)
    # Do not plot axis
    plt.axis("off")
    # Set title
    plt.title("Cell Masks", fontsize=20, color='white')
    # Save figure
    plt.savefig(output_dir+'/cell_segmentation.tif',transparent=True)  
    # Close Fig
    plt.close()
    print("2. Exported cell masks on cell + nuclear signal!", flush=True)
    
    # 3. Cytoplasmic masks on protein signal
    # Create Figure Object
    plt.figure(figsize=(6,6))
    # Plot Raw nuclei with selected percentiles 
    plt.imshow(protein_raw, vmin=np.percentile(protein_raw,vmin_prot), vmax=np.percentile(protein_raw,vmax_prot), cmap='gray')
    # Plot Nuclear segmentation masks
    plt.imshow(mycmap(masks_cyto_only), alpha=0.5)
    # Do not plot axis
    plt.axis("off")
    # Set title
    plt.title("Cyplasmic Masks on Protein Signal", fontsize=20, color='white')
    # Save figure
    plt.savefig(output_dir+'/cyto_mask_on_protein.tif',transparent=True)
    # Close Fig
    plt.close()
    print("3. Exported cytoplasmic mask on protein signal!", flush=True)
    
    # 4. Labelled cell masks
    # Create figure object
    plt.figure(figsize=(6,6))
    # Plot cell masks
    plt.imshow(mycmap(masks_cell_fil))
    # Loop over cell masks
    for i in range(0,len(props_protein_in_cell)):
    
        # Get coordinates of cytoplasm centroids
        coord=props_protein_in_cell[i].centroid
    
        # Annotate the id of the cells
        plt.annotate(props_protein_in_cell[i].label,(coord[1],coord[0]),fontsize=3,color='white')
    # Set title
    plt.title("Labelled Cell Masks", fontsize=20, color='white')
    # Do not plot figure axis
    plt.axis('off')
    # Save figure
    plt.savefig(output_dir+'/cell_segmentation_label.tif',transparent=True)
    # Close Fig
    plt.close()
    print("4. Exported labelled cell masks!", flush=True)
    
    # 5. Labelled nuclear masks
    # Create figure object
    plt.figure(figsize=(6,6))
    # Plot cell masks
    plt.imshow(mycmap(masks_nu_in_cell))
    # Loop over cell masks
    for i in range(0,len(props_protein_in_nu)):
    
        # Get coordinates of cytoplasm centroids
        coord=props_protein_in_nu[i].centroid
    
        # Annotate the id of the cells
        plt.annotate(props_protein_in_nu[i].label,(coord[1],coord[0]),fontsize=3,color='white')
    
    # Set title
    plt.title("Labelled Nuclear Masks", fontsize=20, color='white')
    # Do not plot figure axis
    plt.axis('off')
    # Save figure
    plt.savefig(output_dir+'/nuclear_segmentation_label.tif',transparent=True)
    # Close Fig
    plt.close()
    print("5. Exported labelled nuclear masks!", flush=True)
    
    # 6. Table of properties
    header = ['#id', 'centroid_nu_x', 'centroid_nu_y', 'area_cell', 'area_nu', 'area_cyto', 'mean_int_prot_cell', 'mean_int_prot_nu', 'mean_int_prot_cyto', 'cell_major_axis', 'cell_minor_axis', 'cell_aspect_ratio', 'nucleus_major_axis', 'nucleus_minor_axis', 'nucleus_aspect_ratio']
    
    # open the file in the write mode
    with open(output_dir+'/properties.csv', 'w') as f:
        # create the csv writer
        writer = csv.writer(f)
        # write a row to the csv file
        writer.writerow(header)
        for i in range(0,len(props_protein_in_nu)):
                if props_protein_in_cell[i].axis_minor_length == 0 and props_protein_in_nu[i].axis_minor_length == 0:
                    writer.writerow([props_protein_in_nu[i].label,props_protein_in_nu[i].centroid[0],props_protein_in_nu[i].centroid[1],props_protein_in_cell[i].area*pixel_size*pixel_size,props_protein_in_nu[i].area*pixel_size*pixel_size,
                                props_protein_in_cyto[i].area*pixel_size*pixel_size,props_protein_in_cell[i].intensity_mean,props_protein_in_nu[i].intensity_mean,props_protein_in_cyto[i].intensity_mean,
                                props_protein_in_cell[i].axis_major_length,props_protein_in_cell[i].axis_minor_length,"NA",
                                props_protein_in_nu[i].axis_major_length,props_protein_in_nu[i].axis_minor_length,"NA"])
                elif props_protein_in_cell[i].axis_minor_length == 0:
                    writer.writerow([props_protein_in_nu[i].label,props_protein_in_nu[i].centroid[0],props_protein_in_nu[i].centroid[1],props_protein_in_cell[i].area*pixel_size*pixel_size,props_protein_in_nu[i].area*pixel_size*pixel_size,
                                props_protein_in_cyto[i].area*pixel_size*pixel_size,props_protein_in_cell[i].intensity_mean,props_protein_in_nu[i].intensity_mean,props_protein_in_cyto[i].intensity_mean,
                                props_protein_in_cell[i].axis_major_length,props_protein_in_cell[i].axis_minor_length,"NA",
                                props_protein_in_nu[i].axis_major_length,props_protein_in_nu[i].axis_minor_length,(props_protein_in_nu[i].axis_major_length/props_protein_in_nu[i].axis_minor_length)])
                elif props_protein_in_nu[i].axis_minor_length == 0:
                    writer.writerow([props_protein_in_nu[i].label,props_protein_in_nu[i].centroid[0],props_protein_in_nu[i].centroid[1],props_protein_in_cell[i].area*pixel_size*pixel_size,props_protein_in_nu[i].area*pixel_size*pixel_size,
                                props_protein_in_cyto[i].area*pixel_size*pixel_size,props_protein_in_cell[i].intensity_mean,props_protein_in_nu[i].intensity_mean,props_protein_in_cyto[i].intensity_mean,
                                props_protein_in_cell[i].axis_major_length,props_protein_in_cell[i].axis_minor_length,(props_protein_in_cell[i].axis_major_length/props_protein_in_cell[i].axis_minor_length),
                                props_protein_in_nu[i].axis_major_length,props_protein_in_nu[i].axis_minor_length,"NA"])
                else:
                    writer.writerow([props_protein_in_nu[i].label,props_protein_in_nu[i].centroid[0],props_protein_in_nu[i].centroid[1],props_protein_in_cell[i].area*pixel_size*pixel_size,props_protein_in_nu[i].area*pixel_size*pixel_size,
                                props_protein_in_cyto[i].area*pixel_size*pixel_size,props_protein_in_cell[i].intensity_mean,props_protein_in_nu[i].intensity_mean,props_protein_in_cyto[i].intensity_mean,
                                props_protein_in_cell[i].axis_major_length,props_protein_in_cell[i].axis_minor_length,(props_protein_in_cell[i].axis_major_length/props_protein_in_cell[i].axis_minor_length),
                                props_protein_in_nu[i].axis_major_length,props_protein_in_nu[i].axis_minor_length,(props_protein_in_nu[i].axis_major_length/props_protein_in_nu[i].axis_minor_length)])
    
    # Close file
    f.close();
    print("6. Exported table of properties!", flush=True)
    
