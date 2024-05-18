as_regulon_network = function(regulons){

    regulators = regulons[['regulator']] %>% unique()
    regulons = sapply(regulators, function(regulator_oi){
            X = regulons %>%
                filter(regulator %in% regulator_oi)
            X = list(
                tfmode = setNames(X[['tfmode']], X[['target']]),
                likelihood = X[['likelihood']]
            )
            return(X)
        }, simplify=FALSE)
    
    return(regulons)
}


load_networks = function(network_path, n_tails="two", patt=NULL){
    if (file.exists(network_path) && !dir.exists(network_path)){
        # network_path is a file, we load only that network (we'll run regular VIPER)
        network_files = list(network_path)
    }else if (dir.exists(network_path)){
        # network_path is a directory, we load all networks contained (we'll run metaVIPER)
        network_files = list.files(network_path, pattern=patt, full.names=TRUE)
    }else {
        stop("Invalid network_path.")
    }
    
    networks = sapply(network_files, function(network_file){
        print(network_file)
        network = read_tsv(network_file)
        if (nrow(network)>1 & n_tails=="one"){
            network = network %>% mutate(tfmode=abs(tfmode))
        }
        network = as_regulon_network(network)
        return(network)
    }, simplify=FALSE)
    
    # drop regulons that cannot be used
    networks = sapply(networks, function(network){
        to_keep = sapply(network, function(x){ length(x[[1]]) }) >= 25 # viper's default minsize
        network = network[to_keep]
        return(network)
    }, simplify=FALSE)
    
    # drop networks that cannot be used
    to_keep = sapply(networks, length)>1
    networks = networks[to_keep]
    
    return(networks)
}