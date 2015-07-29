
function nam = get_name(com,detailed) 

nam = com.name;   

%%%% SGE support
if isdeferred(com), nam = get_name(com.deferred, nam, com.is_data); end
%%%%
