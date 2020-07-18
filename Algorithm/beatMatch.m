function matches = beatMatch(beats1, beats2, Fs)

% Find matching beats in two diffrent signals
%
%   matches = beatMatch(beats1, beats2, Fs)
%
% Inputs:
%   beats1: beat annotations from the first signal (Nx1)
%   beats2: beat annotations from the second signal (Mx1)
%   Fs: sample frequency of the signals
%
% Outputs:
%   matches: a vector Nx1, containing 1 or 0 for every beat. 1 = match, 
%   0 = no match 
% 
% Written by: Linda Eerikainen, 26 March, 2015
% Last modified: tolerance changed 15 July, 2015, by Linda Eerikäinen

n_samples = length(beats1);
tolerance = 0.08*Fs;
matches = zeros(n_samples,1);


for i = 1:n_samples;
   beat = beats1(i);
   match_beat = beats2 >= beat-tolerance & beats2 <= beat+tolerance;
   if nnz(match_beat) == 1
       matches(i) = 1;
   end

end




end

