function synFinal=get_syn_pos_on_axons_with_equal_interdist(pos,dA,synResolution)
% to interpolate equal length along axon
% pos is a Nx3 matrix which contains all axon positons.
axonPosCell = classify_it_into_multi_thread(pos,dA);
synFinal = [];
for num = 1:length(axonPosCell)
    pos = axonPosCell{num,1};
    numPoint = size(pos,1);
    % for i=1:3
    % pos(:,i)=pos(:,i)-pos(1,i)*ones(numPoint,1);% set the first point as origin
    
    interDist = zeros(numPoint-1,1);
    for i=1:numPoint-1
        interDist(i)=norm(pos(i+1,:)-pos(i,:));
    end
    
    numSyn = floor(sum(interDist)/synResolution);
    if numSyn==0
        synPosition = (pos(1,:)+pos(end,:))/2;
    else
        synPosition = zeros(numSyn,3);
        indexStart = 1;
        distSum = 0;
        for i=1:numSyn
            
            while (distSum<synResolution)&&(indexStart<=numPoint)
                distSum = distSum + interDist(indexStart);
                indexStart=indexStart+1;
            end
            remainDist = synResolution+interDist(indexStart-1)-distSum;
            % apply linear interpolation
            ratio = remainDist/interDist(indexStart-1);
            synPosition(i,:) = pos(indexStart-1,:)+ratio*(pos(indexStart,:)-pos(indexStart-1,:));
            distSum=interDist(indexStart-1)-remainDist;
        end
        
    end
    synFinal = cat(1,synFinal,synPosition);
end


end

function posCell = classify_it_into_multi_thread(pos,dA)
% classify axons or dendrites into multiple threads, each thread is connected only themselves.
[row,col] = find(dA);
dif = row - col;
matching = [row,col,dif];
redundantIndex = find(matching(:,3)~=1);
matching(redundantIndex,:)=[];
orderNum = 1:size(matching,1);
orderNum = transpose(orderNum);
matching(:,3) = matching(:,2) - orderNum;
uniNumber = unique(matching(:,3));
posCell = cell(length(uniNumber),1);
for i = 1:length(posCell)
    index = find(matching(:,3) == uniNumber(i));
    posCell{i,1} = pos(matching(index,2),:); % 2 denotes parent nodes
end
end