function makeDir(originalDir)
    % "/"でディレクトリーを分割
    dirs = split(originalDir, "/");
    rootDir = fullfile(dirs{1}, dirs{2});
    
    % 現在のディレクトリをrootDirに設定
    currDir = rootDir;
    
    % それぞれのディレクトリに対して
    for i = 3:length(dirs)
        % 次のディレクトリパスを作成
        currDir = fullfile(currDir, dirs{i});
        
        % ディレクトリが存在するか確認f
        if ~isfolder(currDir)
            % 存在しない場合は新しく作成する
            mkdir(currDir);
        end
    end
end