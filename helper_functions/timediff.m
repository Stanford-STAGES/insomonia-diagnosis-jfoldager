function times = timediff(startdatevec,enddatevec)
    for i = 1:size(startdatevec,1)
        if etime(enddatevec(i,:),startdatevec(i,:)) < 0
           enddatevec(i,3) = enddatevec(i,3) + 1;
        end
    end
    times = etime(enddatevec,startdatevec);
end