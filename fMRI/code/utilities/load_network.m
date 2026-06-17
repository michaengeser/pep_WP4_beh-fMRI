function [net, cfg] = load_network(cfg)
% load network
if strcmp(cfg.dnn,'googlenet_places365')
    net=googlenet('weights','places365');

    % get layers
    if strcmp(cfg.layer_type,'early_mid_late')
        cfg.loi = [8,68,139];
    elseif strcmp(cfg.layer_type,'late')
        cfg.loi= 139;
    elseif strcmp(cfg.layer_type,'early_late')
        cfg.loi= [8,139];
    elseif strcmp(cfg.layer_type,'mid_late')
        cfg.loi= [68,139];
    elseif strcmp(cfg.layer_type,'all_output')
        cfg.loi=[8,25,39,54,68,82,96,110,125,139,142];
    elseif strcmp(cfg.layer_type,'all')
        cfg.loi = 1:144;
    end

elseif contains(cfg.dnn,'vgg')

    if strcmp(cfg.dnn,'vggnet16_places365')
        load(fullfile(pwd,'..', 'vgg16_places365/vggnet16_places365.mat'));

    elseif strcmp(cfg.dnn,'vgg16_imagenet')
        load(fullfile(pwd,'..', 'vgg16_imagenet/vgg16.mat'));
        net=vgg16;
        clear vgg16
    end

    % get layers
    if strcmp(cfg.layer_type,'early_mid_late')
        cfg.loi = [4,21,36];
    elseif strcmp(cfg.layer_type,'late')
        cfg.loi = 36;
    elseif strcmp(cfg.layer_type,'early_late')
        cfg.loi= [4,36];
    elseif strcmp(cfg.layer_type,'mid_late')
        cfg.loi = [21,36];
    elseif strcmp(cfg.layer_type,'all_output')
        %Get the Conv and FC Layers
        lx=0;
        for layer=1:length(net.Layers)
            if strfind(net.Layers(layer).Name,'conv')==1
                lx=lx+1;
                cfg.loi(lx)=layer;
            elseif strfind(net.Layers(layer).Name,'fc')==1
                lx=lx+1;
                cfg.loi(lx)=layer;
            end
        end
    elseif strcmp(cfg.layer_type,'all')
        cfg.loi = 1:41;
    end
end
end