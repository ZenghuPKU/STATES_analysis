function Ripley_k_value = ripley_k_function_edge_corr_v1(pos_points,pos_border, t, param)
%%%  Ripley-k function with border correction
% param.correction = 0 ==> no correction
% param.correction = 1 ==> Correction with renormalization of border  spot count
% param.correction = 2 ==> Correction with deletion of border spot

%% Get some parameters
n_RNA     = param.n_RNA;
mask_cell = param.mask_cell;

%% Some house-keeping
Ripley_k_value = zeros(length(t),1);
lambda         = n_RNA/param.area_cell ;


%% Get the distace matrix and the number of neighbour point
pos          = [pos_points.x  pos_points.y];
dist_mat     = distmat(pos);

% Remove duplicate positions in outlines ... 
border_pos = [];
border_pos(:,1) = pos_border.x;
border_pos(:,2) = pos_border.y;
border_pos = unique(border_pos, 'rows');


%=== Correction with normalized counts at border
if param.correction == 1
   
    % Get the distance from the border
    dist_cell = zeros(size(pos,1),1);
    for i_point = 1:n_RNA
        dist_cell(i_point) = abs(p_poly_dist(pos(i_point,1),pos(i_point,2),border_pos(:,1),border_pos(:,2))) ;
    end
    
    for i_dist = 1:length(t)
        
        dist_mat_ind  = double(dist_mat < t(i_dist) ) ;
        ind_dist      = sum(dist_mat_ind) - 1 ;
        ind_corrected = dist_cell < t(i_dist) ;
        ind_corrected = find(ind_corrected == 1);
        
        %% Do the edge correction
        weigth      = ones( 1,length(dist_cell)) ;
        
        SE = strel('disk',  t(i_dist), 8);
        mask_circle = SE.getnhood();
        
        size_circle = length(find(mask_circle == 1)) ;
        
        
        for i = 1:length(ind_corrected)
            
            i_spot       = ind_corrected(i) ;
            dum_spot     = round([pos_points.y(i_spot)   pos_points.x(i_spot)]);
                    
            coord_crop_x_min = max(1,dum_spot(1) - (size(mask_circle,1) -1)/2);
            coord_crop_x_max = min(dum_spot(1) + (size(mask_circle,1) -1)/2,size(mask_cell,2));
            coord_crop_y_min = max(dum_spot(2) - (size(mask_circle,2) -1)/2,1);
            coord_crop_y_max = min(dum_spot(2) + (size(mask_circle,2) -1)/2,size(mask_cell,1));
            
            coord_crop_x_min_circle = max(1, (size(mask_circle,1) -1)/2 - dum_spot(1)+2);
            coord_crop_x_max_circle = min(size(mask_circle,1), (size(mask_circle,1) -1)/2 + size(mask_cell,2) - dum_spot(1)+1);
            
            coord_crop_y_min_circle = max(1, (size(mask_circle,2) -1)/2 - dum_spot(2)+2);
            coord_crop_y_max_circle = min(size(mask_circle,2),(size(mask_circle,2) -1)/2 + size(mask_cell,1) - dum_spot(2)+1);
            
            crop_mask        = mask_cell(coord_crop_y_min:coord_crop_y_max,coord_crop_x_min:coord_crop_x_max);
            crop_circle      = mask_circle(coord_crop_y_min_circle:coord_crop_y_max_circle,coord_crop_x_min_circle:coord_crop_x_max_circle);
            weigth(i_spot)   = length(find(crop_circle.*crop_mask == 1))/size_circle;
            
        end
        
        
        %==== Calculate ripley-k function
        K                      = sum(ind_dist./weigth)/n_RNA ;
        Ripley_k_value(i_dist) = K/lambda ;
        
    end
    %=== Correction by ignoring spots at the border
elseif param.correction == 2
    
    dist_cell = zeros(size(pos,1),1);
    for i_point = 1:n_RNA
        dist_cell(i_point) = abs(p_poly_dist(pos(i_point,1),pos(i_point,2), pos_border.x,pos_border.y)) ;
    end
    
    for i_dist = 1:length(t)
        
        dist_mat_ind  = double(dist_mat < t(i_dist) ) ;
        ind_dist      = sum(dist_mat_ind) - 1 ;
        ind_spot      = dist_cell > t(i_dist);
        ind_spot      = find(ind_spot == 1);
        
        value_spot    = ind_dist(ind_spot) ;
        
        K                      = sum(value_spot)/length(ind_spot) ;
        Ripley_k_value(i_dist) = K/lambda ;
        
        
    end
    
    
%=== No correction    
elseif param.correction == 0
    
    for i_dist = 1:length(t)
        
        dist_mat_ind  = double(dist_mat < t(i_dist) ) ;
        ind_dist      = sum(dist_mat_ind) - 1 ;
        K                      = sum(ind_dist)/n_RNA ;
        Ripley_k_value(i_dist) = K/lambda ;
        
    end  
end


