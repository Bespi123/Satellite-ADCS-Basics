% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Copyright 2016 Mayuresh Sarpotdar
% % 
% %    Licensed under the Apache License, Version 2.0 (the "License");
% %    you may not use this file except in compliance with the License.
% %    You may obtain a copy of the License at
% % 
% %        http://www.apache.org/licenses/LICENSE-2.0
% % 
% %    Unless required by applicable law or agreed to in writing, software
% %    distributed under the License is distributed on an "AS IS" BASIS,
% %    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% %    See the License for the specific language governing permissions and
% %    limitations under the License.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%Modified by bespi123: 26/10/2023

function [seg, real_id] = select_region_fov(star_cat, fov, u, R)
u = u/norm(u);
ux_cent = u(1); uy_cent = u(2); uz_cent = u(3);

count = 0;
for i=1:length(star_cat)
    d = acosd( ux_cent*star_cat(i,4) + uy_cent*star_cat(i,5) + uz_cent*star_cat(i,6));
    if (d<fov/2)
        count = count + 1;
        seg(count,1) = count;
        seg(count,2:3) = star_cat(i,2:3);
        seg(count,4:6) = star_cat(i,4:6); 
        seg(count,7:9) = (R*star_cat(i,4:6)')'; % unit vectors in the rotated coordinate system
%         ra_star_seg = star_cat(i,2) - ra;
%         dec_star_seg = star_cat(i,3) - dec;
%         seg(count,4) = cosd(dec_star_seg)*cosd(ra_star_seg);
%         seg(count,5) = sind(dec_star_seg); % cosd(dec_star_seg)*sind(ra_star_seg)
%         seg(count,6) = cosd(dec_star_seg)*sind(ra_star_seg); % sind(dec_star_seg);
        real_id(count,1) = i;
    end
end

if(count==0)
    seg = zeros(1,6);
    real_id = 0;
end
