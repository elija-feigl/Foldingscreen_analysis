function [parsed_data, warnings] = parse_gel_info(filepath, log_file)
% @ step1 via compute_profiles
% parse the gel_info.txt file

    %% read line by line
    warnings = false;
    disp(filepath)
    fileID = fopen(filepath);
    tmp = textscan(fileID,'%s','CommentStyle','#', 'Delimiter', '\n');
    fclose(fileID);
    lines = tmp{1};
    parsed_data.filepath = filepath;
    parsed_data.log_file = log_file;
    
    logfile_ID = fopen(log_file,'a');
    fprintf(logfile_ID,'%s\n', ['Parsing file: ' filepath]);
    disp(['Parsing file: ' filepath])
    
    
    info = ["user", "project", "design_name", "date", "scaffold_type", "lattice_type", "tem_verified" ,"comment", "scaffold_concentration", "staple_concentration","comment", "Lane_" "Gelsize", "Agarose_concentration","Staining", "Mg_concentration" ,"Voltage", "Running_time","Cooling"];
    for i=1:length(lines)
        if ~isempty(lines{i}) % if line containes text, e.i. is not empty
            seg = split(lines{i}, {'=', ':'}); % split at = or :
            if length(seg)==1 % there was no = or : in the line
                disp(['Warning. Random line detected. It will be ignored. Check gel_info. Line: ' lines{i}])
                fprintf(logfile_ID,'%s\n', ['Warning. Random line detected. It will be ignored. Check gel_info. Line: ' lines{i}]);

                warnings = true;
            else
                index_comment = strfind(seg{2}, '#');
                if ~isempty(index_comment)
                    seg{2} = seg{2}(1:index_comment(1)-1); % remove all characters afer comment 
                end
                seg{2} = strtrim(seg{2});
            end
            
            
            key = strrep(seg{1}, ' ', '');
            
            if contains(key, "Lane_")
                l = str2num(erase(key, "Lane_"));
                key = "Lane_";
            end
            if  any (strcmpi (info, key))
                if any (strcmpi(["scaffold_concentration", "staple_concentration"], key))
                    parsed_data.(lower(key)) = str2num(strtrim(seg{2}));
                elseif key == "Lane_"
                    parsed_data.lanes{l} = strtrim(seg{2});
                else
                    parsed_data.(lower(key)) = strtrim(seg{2});
                end
            end  
        end
    end
    %disp(['File ' filepath ' parsed.'])
    
    %% go through lanes
    
    advanced = false;
    for i=1:length(parsed_data.lanes)
            %disp(parsed_data.lanes{i})
            if strcmp(parsed_data.lanes{i}(1), '{')
                advanced = true;
            end
    end
    
    if advanced
        parsed_data.lanes_unparsed = parsed_data.lanes; % save un-parsed data
        for i=1:length(parsed_data.lanes)
            %disp(parsed_data.lanes{i})
            if strcmp(parsed_data.lanes{i}(1), '{')
                tmp = strsplit(parsed_data.lanes{i}(2:end-1), ',');
                parsed_data.lanes{i} = tmp{1};
            end
        end
    end
    
    %% put each lane into its catagory: (ladder, scaffold, mono)
    parsed_data.lanes = string(parsed_data.lanes);
    indices_list = 1:length(parsed_data.lanes); 
    ladder_indices = find(contains(parsed_data.lanes, "ladder"));
    scaff_indices = find(contains(parsed_data.lanes, "scaff"));
    
    if ladder_indices
        parsed_data.species.ladder.type = "ladder";
        parsed_data.species.ladder.indices = ladder_indices ;
        parsed_data.species.ladder.lane_name = parsed_data.lanes(ladder_indices);

        
    else
        parsed_data.species.ladder = [];
        parsed_data.species.ladder.type = "none";
        parsed_data.species.ladder.indices = [];
    end 
    
    if scaff_indices
        parsed_data.species.scaffold.type = "scaffold";
        parsed_data.species.scaffold.indices = scaff_indices ;
        parsed_data.species.scaffold.lane_name = parsed_data.lanes(scaff_indices);
    else
        parsed_data.species.scaffold = [];
        parsed_data.species.scaffold.type = "none";
        parsed_data.species.scaffold.indices = [];
    end
    
    indices_list([ladder_indices scaff_indices]) = [];
    parsed_data.species.mono.type = "mono";
    
    
    parsed_data.species.mono.indices = indices_list;
    parsed_data.species.mono.lane_name = parsed_data.lanes(indices_list);
    
    species = parsed_data.species;
    parsed_data.loop_indices = cat(2, species.ladder.indices , species.scaffold.indices , species.mono.indices);
    
    
    
   fclose(logfile_ID);

end








