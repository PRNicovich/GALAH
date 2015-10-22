%%%%%%%%%%%%%%%%%%%%
% Clean up sk in Calc_some_ROIs

    function skClean = CleanUpSk(sk)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
                    
        threshold_val = get(handles.handles.ROI_thresh_slide_hand, 'Value');
        frame_here = str2double(get(handles.handles.ROI_frame_slide_box, 'String'));
        image_here = handles.Img_stack(:,:, frame_here, handles.Primary_channel);
        
        
        if handles.Primary_channel == 1;
            min_max_here = handles.Min_max_left;
            display_here = handles.Display_range_left;
        elseif handles.Primary_channel == 2;
            min_max_here = handles.Min_max_right;
            display_here = handles.Display_range_right;
        end
        
        min_max_here = [min(min_max_here(:,1)), max(min_max_here(:,2))];
        
        scale_frame = (image_here - min_max_here(1))./(min_max_here(2) - min_max_here(1));

        [sizx, sizy] = size(image_here);

        
        branch_im = bwmorph(sk,'branchpoints',1);
        cycles_im = imfill(sk, 4, 'holes') - sk;
        
        bpoints = sum(branch_im(:));
        cpoints = sum(cycles_im(:));
        
        while (bpoints + cpoints) > 0
            
            if bpoints > 0
                
                % For each bpoints...
                bpoints_list = find(branch_im);
                
                bpdil = bwmorph(branch_im, 'dilate');
                bp = sk;
                bp(bpdil) = 0;
                
                CC = bwconncomp(bp);
                segLabels = labelmatrix(CC);
                
                for k = 1:numel(bpoints_list);
                    
                    % Pull out bpoints(k)
                    phere = bpoints_list(k);
                    [psubx, psuby] = ind2sub([sizx sizy], phere);
                    
                    % Find which segments in bp connect to bpoint(k)
                    conn8 = segLabels((psubx + [-2:2]), (psuby + [-2:2]));
%                     segConns = unique(conn8(conn8 > 0));
                    segAll = unique(conn8(conn8>0));
                    
                    % Calc mean (total?) brightness of connected line segments
                    segMeans = zeros(numel(segAll), 1);
                    segTotal = zeros(numel(segAll), 1);
                    segNumel = zeros(numel(segAll), 1);
                    
                    for m = 1:numel(segAll)
                        
                        segMeans(m) = mean(image_here(CC.PixelIdxList{segAll(m)}));
                        segTotal(m) = sum(image_here(CC.PixelIdxList{segAll(m)}));
                        segNumel(m) = numel(image_here(CC.PixelIdxList{segAll(m)}));
                        
                    end
                    
                    % Zero out all but the top two segments in that branch point in sk
                    sortSegs = sort((segMeans), 'descend');
                    toClear = segAll(ismember(segMeans, sortSegs(3:end)));
                    
                    for m = 1:numel(toClear)
                        
                        sk(CC.PixelIdxList{toClear(m)}) = 0;
                        
                    end
                    
                end
                
                sk = bwmorph(sk, 'thin', Inf);
                sk = bwmorph(sk, 'spur');
                sk = bwmorph(sk, 'clean');
            end
            
            if cpoints > 0
                
                cycles_im = imfill(sk, 4, 'holes') - sk;
                
                bdil = bwmorph(cycles_im, 'dilate');
                cbranch = bwmorph(sk, 'branchpoints');
                sk(bdil & ~cbranch) = 0;
                sk = bwmorph(sk, 'spur');
                sk = bwmorph(sk, 'clean');
                
            end
            
            bpoints = sum(sum(bwmorph(sk,'branchpoints',1)));
            cpoints = sum(sum(imfill(sk, 4, 'holes') - sk));
        end
        
        
        skClean = sk;
    end