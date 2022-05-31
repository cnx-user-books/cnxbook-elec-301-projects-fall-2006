function SECBIGGER = framedifference (SECTIONS)

ssize = size(SECTIONS);

for ii = ssize(1):-1:2
    FREQDIFF(ii-1,:) = SECTIONS(ii-1,:)-SECTIONS(ii,:);
end

for jj = 2:ssize(2)
    SECDIFF(:,jj-1) = FREQDIFF(:,jj) - FREQDIFF(:,jj-1);
end

for ii = 1:ssize(1)-1
    for jj = 1:ssize(2)-1
        if SECDIFF(ii,jj) > 0
            SECBIGGER(ii,jj) = 1;
        end
    end
end
