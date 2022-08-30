function [runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG(whoisrunning,print_info)
%add_paths_matlab(print_info)
%  
if(nargin == 0)
    print_info = 1;
    [~,whoisrunning] = system('whoami')
elseif(nargin == 1)
    print_info = 1;
end

switch strtrim(whoisrunning)
    case 'duip66033\lpzmi_local'
        main_path                       = 'D:/Dropbox';
        runpath                         = fullfile(main_path,'fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/');
        code_path.fieldtrippath         = fullfile(main_path,'big_files_shared','fieldtrip-20170827');
        %code_path.fieldtrippath=fullfile('D:\toolbox','fieldtrip-20180711');
        code_path.my_functions          = fullfile(runpath,'my_functions');
        %session_path.raw=fullfile('H:\MEGEYE\CTF_DATA');
        
        session_path.raw        = fullfile(main_path,'big_files_shared/CTF_DATA');
        session_path.rawet      = fullfile(main_path,'big_files_shared/ET_DATA');
        session_path.matfiles   = fullfile(main_path,'big_files_shared/MAT_DATA');
        session_path.out        = fullfile(main_path,'OUT\Analysis_ET_Beh');
        
    case 'ad\lpajg1' % Joac Notts
        main_path                       = 'C:\Users\lpajg1\OneDrive - The University of Nottingham\MEGEYEHS\';
        runpath                         = fullfile(main_path, 'Matlab_analysis\');
        code_path.fieldtrippath         = fullfile(main_path, 'Matlab_analysis\fieldtrip-20220104\');
        code_path.mrTools               = fullfile(main_path, 'Matlab_analysis\mrTools-master\');
        code_path.my_functions          = fullfile(runpath,'my_functions');
        local_path_beh                  = fullfile(main_path,'Psychopy_experiment\'); 
        
        session_path.meg        = fullfile(main_path,'DATA\CTF_DATA');
        session_path.et      = fullfile(main_path,'DATA\ET_DATA');
        session_path.matfiles   = fullfile(main_path,'DATA\MAT_DATA');
        session_path.behav = fullfile(main_path,'DATA\BH_DATA');
        session_path.mri        = fullfile(main_path,'DATA\MRI_DATA');
        session_path.opt        = fullfile(main_path,'DATA\OPT_DATA');
        session_path.out        = fullfile(main_path,'OUT\Analysis_ET_Beh');
        session_path.positions_file  = fullfile(main_path,'DATA','pos_items_210.mat'); %IMPORTANT: This file...
        %contains the annotated positions of all the 210 stimuli used.
        %Currently it only has information on 20-30 of them
        session_path.images_folder   = local_path_beh;
    
    case 'laptop-5i5qsv76\joaco' % Joac Asus
        main_path                       = 'C:\Users\joaco\OneDrive - The University of Nottingham\MEGEYEHS\';
        runpath                         = fullfile(main_path, 'Matlab_analysis\');
        code_path.fieldtrippath         = fullfile(main_path, 'Matlab_analysis\fieldtrip-20220104\');
        code_path.mrTools               = fullfile(main_path, 'Matlab_analysis\mrTools-master\');
        code_path.my_functions          = fullfile(runpath,'my_functions');
        local_path_beh                  = fullfile(main_path,'Psychopy_experiment\'); 
        
        session_path.meg        = fullfile(main_path,'DATA\CTF_DATA');
        session_path.et      = fullfile(main_path,'DATA\ET_DATA');
        session_path.matfiles   = fullfile(main_path,'DATA\MAT_DATA');
        session_path.behav = fullfile(main_path,'DATA\BH_DATA');
        session_path.mri        = fullfile(main_path,'DATA\MRI_DATA');
        session_path.opt        = fullfile(main_path,'DATA\OPT_DATA');
        session_path.out        = fullfile(main_path,'OUT\Analysis_ET_Beh');
        session_path.positions_file  = fullfile(main_path,'DATA','pos_items_210.mat'); %IMPORTANT: This file...
        %contains the annotated positions of all the 210 stimuli used.
        %Currently it only has information on 20-30 of them
        session_path.images_folder   = local_path_beh;

%     case 'root' %'juank-Desktop'
%         main_path                       = '~/Dropbox';
%         runpath                         = fullfile(main_path,'fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/');
%         code_path.fieldtrippath         = fullfile('/home/juank/toolbox/fieldtrip-20170827');
%         code_path.my_functions          = fullfile(runpath,'my_functions');
%         
%         session_path.raw        = fullfile(main_path,'big_files_shared/CTF_DATA');
%         session_path.rawet      = fullfile(main_path,'big_files_shared/ET_DATA');
%         session_path.matfiles   = fullfile(main_path,'big_files_shared/MAT_DATA');
%         session_path.out        = fullfile(main_path,'OUT\Analysis_ET_Beh');

    case 'joaco'
        main_path                       = '/mnt/6a6fd40a-e256-4844-8004-0e60d95969e8/MEGEYEHS/';
        runpath                         = fullfile(main_path,'Matlab/');
        code_path.fieldtrippath         = fullfile('/home/usuario/fieldtrip-20220827');
        code_path.my_functions          = fullfile(runpath,'my_functions');
%         code_path.mrTools               = fullfile(main_path, 'Matlab\mrTools-master\');
%         local_path_beh                  = fullfile(main_path,'Psychopy_experiment\'); 
        
        session_path.meg        = fullfile(main_path,'DATA/CTF_DATA');
        session_path.et         = fullfile(main_path,'DATA/ET_DATA');
        session_path.matfiles   = fullfile(main_path,'DATA/MAT_DATA');
        session_path.behav      = fullfile(main_path,'DATA/BH_DATA');
        session_path.mri        = fullfile(main_path,'DATA/MRI_DATA');
        session_path.opt        = fullfile(main_path,'DATA/OPT_DATA');
        session_path.out        = fullfile(main_path,'OUT/');
        session_path.positions_file  = fullfile(main_path,'DATA','pos_items_210.mat'); %IMPORTANT: This file...
        %contains the annotated positions of all the 210 stimuli used.
        %Currently it only has information on 20-30 of them
%         session_path.images_folder   = local_path_beh;

        fprintf('WARNING: Define path to psychopy experiment to access the images\n')
        
    otherwise
        warning(['User unknown: ' whoisrunning ' , check paths']);
end

% session_path.directories = {''; ''};
session_path.subjname = {'15909001','15912001','15910001','15950001','15911001','11535009'};
session_path.subjcode = {'15909001','15912001','15910001','15950001','15911001','11535009'};
session_path.sessionfilenames = {
    {'15909001\15909001_MatiasIson_20220530_02.ds',  '15909001\15909001_MatiasIson_20220530_03.ds',...
    '15909001\15909001_MatiasIson_20220530_04.ds'};...
    {'15912001\15912001_MatiasIson_20220531_01.ds',  '15912001\15912001_MatiasIson_20220531_02.ds'};...
    {'15910001\15910001_MatiasIson_20220609_01.ds', '15910001\15910001_MatiasIson_20220609_02.ds',...
    '15910001\15910001_MatiasIson_20220609_03.ds', '15910001\15910001_MatiasIson_20220609_04.ds'};...
    {'15950001\15950001_MatiasIson_20220610_01.ds',  '15950001\15950001_MatiasIson_20220610_02.ds',...
    '15950001\15950001_MatiasIson_20220610_03.ds',  '15950001\15950001_MatiasIson_20220610_04.ds',...
    '15950001\15950001_MatiasIson_20220610_05.ds',  '15950001\15950001_MatiasIson_20220610_06.ds'};...
    {'15911001\15911001_MatiasIson_20220610_01.ds'};...
    {'11535009\11535009_MatiasIson_20220613_01.ds'}
    };
session_path.behavfilenames = {
    '15909001\15909001_hybrid_search_builder_code_2022_May_30_1628.csv';...
    '15912001\15912001_hybrid_search_builder_code_2022_May_31_1009.csv';...
    '15910001\15910001_hybrid_search_builder_code_2022_Jun_09_1507.csv';...
    '15950001\15950001_hybrid_search_builder_code_2022_Jun_10_1423.csv';...
    '15911001\15911001_hybrid_search_builder_code_2022_Jun_10_1618.csv';...
    '11535009\11535009_hybrid_search_builder_code_2022_Jun_13_1035.csv'
    };   

session_path.local_path_edf = {
    '15909001\et_data_15909001.edf';...
    '15912001\et_data_15912001.edf';...
    '15910001\et_data_15910001.edf';...
    '15950001\et_data_15950001.edf';...
    '15911001\et_data_15911001.edf';...
    '11535009\et_data_11535009.edf'
    }; 

if(print_info)
    disp(whoisrunning)
    disp(code_path)
    disp(runpath)
end

end

