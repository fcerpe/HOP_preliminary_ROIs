function mni = roi_getMNIcoords(subs)

% for each subject in subs return the coords

% VWFA = -44 -56 -16, keyword: visual words
%    lLO  = -47 -70 -5
%    rLO  =  47 -70 -5, keyword: objects

% mni{sub_number}(contrast number, coordinates)
% e.g. mni{1}= [-46 64 4]
mni = cell(length(subs),1);

for i = 1:length(subs)

    switch subs{i}
        
        case '001' 
            mni{i}(1, 1:3) = [-46.80 -70.40 -12.80]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-46.80 -67.80  00.20]; %  LO LEFT
            mni{i}(4, 1:3) = [ 49.40 -62.80 -05.00]; %  LO RIGHT

        case '002'
            mni{i}(1, 1:3) = [-49.40 -57.40 -12.80]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-46.80 -62.60 -07.60]; % LO LEFT
            mni{i}(4, 1:3) = [ 46.80 -70.40 -05.00]; % LO RIGHT

        case '003'
            mni{i}(1, 1:3) = [-49.40 -60.00 -20.60]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-39.00 -73.00 -02.40]; % LO LEFT
            mni{i}(4, 1:3) = [ 44.20 -80.80 -05.00]; % LO RIGHT

        % Control for which localizer did not work
        case '004'
            mni{i}(1, 1:3) = [NaN NaN NaN];  % VWFA FRE
            mni{i}(2, 1:3) = [NaN NaN NaN];  % VWFA BRA
            mni{i}(3, 1:3) = [NaN NaN NaN];  % LO LEFT
            mni{i}(5, 1:3) = [NaN NaN NaN];  % LO RIGHT

        case '005' 
            mni{i}(1, 1:3) = [-65.00 -41.80  05.40]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-44.20 -78.20 -07.60]; % LO LEFT
            mni{i}(4, 1:3) = [ 36.40 -83.40 -10.20]; % LO RIGHT

        case '006'
            mni{i}(1, 1:3) = [-46.80 -60.00 -12.80]; % VWFA FRE
            mni{i}(2, 1:3) = [-49.40 -62.60 -12.80]; % VWFA BRA
            mni{i}(3, 1:3) = [-52.00 -65.20 -00.20]; % LO LEFT
            mni{i}(4, 1:3) = [ 39.00 -67.80 -10.20]; % LO RIGHT

        case '007'
            mni{i}(1, 1:3) = [-52.00 -49.60 -18.00]; % VWFA FRE
            mni{i}(2, 1:3) = [-52.00 -49.60 -18.00]; % VWFA BRA
            mni{i}(3, 1:3) = [-44.20 -65.20 -12.80]; % LO LEFT
            mni{i}(4, 1:3) = [ 46.80 -78.20  00.20]; % LO RIGHT
        
        % Expert subject, no VWFA defined by Braille was found    
        case '008'
            mni{i}(1, 1:3) = [-52.00 -54.80 -18.00]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-41.60 -78.20 -02.40]; % LO LEFT
            mni{i}(4, 1:3) = [ 44.20 -78.20 -07.60]; % LO RIGHT

        case '009'
            mni{i}(1, 1:3) = [-44.20 -62.60 -12.80]; % VWFA FRE
            mni{i}(2, 1:3) = [-41.60 -62.60 -10.20]; % VWFA BRA
            mni{i}(3, 1:3) = [-49.40 -78.20 -02.40]; % LO LEFT
            mni{i}(4, 1:3) = [ 36.40 -70.40 -12.80]; % LO RIGHT

        case '010' 
            mni{i}(1, 1:3) = [-44.20 -62.60 -12.80]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-49.40 -75.60  00.20]; % LO LEFT
            mni{i}(4, 1:3) = [ 49.40 -65.20 -05.00]; % LO RIGHT

        case '011' 
            mni{i}(1, 1:3) = [-49.40 -60.00 -10.20]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-39.00 -73.00 -10.20]; % LO LEFT
            mni{i}(4, 1:3) = [ 44.20 -65.20 -10.20]; % LO RIGHT

        case '012'
            mni{i}(1, 1:3) = [-44.20 -60.00 -15.40]; % VWFA FRE
            mni{i}(2, 1:3) = [-44.20 -57.40 -18.00]; % VWFA BRA
            mni{i}(3, 1:3) = [-49.40 -75.60 -05.00]; % LO LEFT
            mni{i}(4, 1:3) = [ 49.40 -73.00 -10.20]; % LO RIGHT

        case '013'
            mni{i}(1, 1:3) = [-44.20 -60.00 -12.80]; % VWFA FRE
            mni{i}(2, 1:3) = [-46.80 -62.60 -15.40]; % VWFA BRA
            mni{i}(3, 1:3) = [-49.40 -78.20 -02.40]; % LO LEFT
            mni{i}(4, 1:3) = [ 49.40 -67.80 -12.80]; % LO RIGHT

        % Subject dropped out and did not complete the experiment
        case '015' 
            mni{i}(1, 1:3) = [NaN NaN NaN];  % VWFA FRE
            mni{i}(2, 1:3) = [NaN NaN NaN];  % VWFA BRA
            mni{i}(3, 1:3) = [NaN NaN NaN];  % LO LEFT
            mni{i}(4, 1:3) = [NaN NaN NaN];  % LO RIGHT
        
        case '017' % Expert for which localizer did not work
            mni{i}(1, 1:3) = [NaN NaN NaN];  % VWFA FRE
            mni{i}(2, 1:3) = [NaN NaN NaN];  % VWFA BRA
            mni{i}(3, 1:3) = [NaN NaN NaN];  % LO LEFT
            mni{i}(4, 1:3) = [NaN NaN NaN];  % LO RIGHT

        case '018' 
            mni{i}(1, 1:3) = [-46.80 -54.80 -07.60]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-46.80 -73.00 -05.00]; % LO LEFT
            mni{i}(4, 1:3) = [ 44.20 -70.40 -10.20]; % LO RIGHT

        case '019' 
            mni{i}(1, 1:3) = [-46.80 -47.00 -20.60]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-41.60 -78.20 -05.00]; % LO LEFT
            mni{i}(4, 1:3) = [ 44.20 -73.00 -07.60]; % LO RIGHT

        case '020' 
            mni{i}(1, 1:3) = [-44.20 -57.40 -20.60]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-49.40 -78.20 -05.00]; % LO LEFT
            mni{i}(4, 1:3) = [ 44.20 -65.20 -12.80]; % LO RIGHT

        case '021' 
            mni{i}(1, 1:3) = [-46.80 -49.60 -23.20]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-46.80 -65.20 -12.80]; % LO LEFT
            mni{i}(4, 1:3) = [ 49.40 -67.80 -07.60]; % LO RIGHT

        case '022' 
            mni{i}(1, 1:3) = [-47.70 -64.90 -18.00]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-42.50 -75.30 -10.20]; % LO LEFT
            mni{i}(4, 1:3) = [ 49.20 -67.50 -10.20]; % LO RIGHT

        case '023' 
            mni{i}(1, 1:3) = [-46.80 -49.60 -23.20]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [   NaN    NaN    NaN]; % LO LEFT: Originally -46.80 -78.20  00.20, but cluster is not sufficient
            mni{i}(4, 1:3) = [ 39.00 -67.80 -05.00]; % LO RIGHT
    
        case '024' % Control subject, Drawings v. scrambled contrast did not show LOC
            mni{i}(1, 1:3) = [-41.60 -52.20 -10.20]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [   NaN    NaN    NaN]; % LO LEFT
            mni{i}(4, 1:3) = [   NaN    NaN    NaN]; % LO RIGHT

        case '025' % Control subject, niether [fw] - [sfw] nor [ld] - [sld]
            mni{i}(1, 1:3) = [   NaN    NaN    NaN]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [   NaN    NaN    NaN]; % LO LEFT
            mni{i}(4, 1:3) = [   NaN    NaN    NaN]; % LO RIGHT

        case '026' 
            mni{i}(1, 1:3) = [-52.00 -57.40 -15.40]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-49.40 -73.00 -10.20]; % LO LEFT
            mni{i}(4, 1:3) = [ 44.20 -67.80 -12.80]; % LO RIGHT

        case '027' 
            mni{i}(1, 1:3) = [-54.60 -57.40 -15.40]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-52.00 -75.60  00.20]; % LO LEFT
            mni{i}(4, 1:3) = [ 39.00 -75.60 -05.00]; % LO RIGHT

        case '028' 
            mni{i}(1, 1:3) = [-54.60 -49.60 -12.80]; % VWFA FRE
            mni{i}(2, 1:3) = [   NaN    NaN    NaN]; % VWFA BRA
            mni{i}(3, 1:3) = [-46.80 -73.00 -12.80]; % LO LEFT
            mni{i}(4, 1:3) = [ 39.00 -73.00 -10.20]; % LO RIGHT

    end
end

end