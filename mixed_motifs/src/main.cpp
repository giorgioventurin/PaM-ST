#include "mixed.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_map>

namespace fs = std::filesystem;
using namespace mixed;
namespace {
struct Options {
    Config cfg;
    std::string input, output, templates, null_model = "blocks";
    std::string x_column = "X_centroid", y_column = "Y_centroid";
    std::string type_column = "Cell_Type";
    double block_size = 1000;
    int candidate_pool = 100;
    bool overwrite = false;
};
std::string trim(const std::string& s) {
    auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
std::string lower(std::string s) {
    for (char& c : s) if (c >= 'A' && c <= 'Z') c += 'a' - 'A';
    return s;
}
double number(const std::string& s, const std::string& name) {
    std::size_t used = 0;
    double value;
    try { value = std::stod(s, &used); }
    catch (...) { throw std::runtime_error("Invalid number for " + name + ": " + s); }
    if (used != s.size() || !std::isfinite(value))
        throw std::runtime_error("Expected finite number for " + name + ": " + s);
    return value;
}
int integer(const std::string& s, const std::string& name) {
    double value = number(s, name);
    if (value != std::floor(value) || value < -2147483647.0 || value > 2147483647.0)
        throw std::runtime_error("Expected integer for " + name);
    return static_cast<int>(value);
}
void help() {
    std::cout << R"(PaM-Mixed: recurring mixed cell-type neighborhoods
Usage: pam_mixed --input cells.csv --output results [options]

  --k 10                        Total cells INCLUDING center
  --max-radius 100              Maximum center-to-neighbor distance, input units
  --min-types 2                 Minimum represented cell types
  --max-dominance 0.7           Largest allowed type fraction in every match
  --max-relabels 1              Half L1 count distance, 0, 1 or 2
  --allow-type-turnover         Allow matches to change represented type set
  --min-support 5               Minimum greedy cell-disjoint occurrences
  --max-motifs 20               Maximum reported discovery motifs
  --candidate-pool 100          Support-screened enrichment pool (0 = all)
  --null blocks                global | stratified | blocks
  --block-size 1000             Coordinate tile size for blocks null
  --region-size 1000            Tile size for descriptive spatial dispersion
  --permutations 199            Number of independently seeded permutations
  --threads 1                   Parallel permutation workers
  --seed 37                    Nonnegative random seed
  --alpha 0.05                 Search/validation significance threshold
  --contact-radius 0           Pair-distance threshold; 0 disables edge filter
  --min-cross-edge-fraction 0   Optional mixed-edge fraction requirement
  --templates FILE             Validate fixed motif_composition.csv templates
  --x-column X_centroid        Coordinate column override
  --y-column Y_centroid
  --type-column Cell_Type
  --overwrite                  Replace generated files in existing output dir

Optional columns: Cell_ID, Sample, Stratum, Component. Neighborhoods never
cross Sample or Component; Stratum restricts shuffles, not neighborhoods.
All nulls preserve Sample x Component cell-type counts; stratified also preserves
Stratum, blocks also preserves coordinate-tile counts. These are random-
label nulls and do not preserve within-stratum homotypic clustering.
See README.md for statistical interpretation and held-out validation.
)";
}
Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") { help(); std::exit(0); }
        if (arg == "--overwrite") { o.overwrite = true; continue; }
        if (arg == "--allow-type-turnover") { o.cfg.same_types = false; continue; }
        if (i + 1 == argc) throw std::runtime_error("Missing value for " + arg);
        const std::string val = argv[++i];
        if (arg == "--input") o.input = val;
        else if (arg == "--output") o.output = val;
        else if (arg == "--templates") o.templates = val;
        else if (arg == "--null") o.null_model = val;
        else if (arg == "--x-column") o.x_column = val;
        else if (arg == "--y-column") o.y_column = val;
        else if (arg == "--type-column") o.type_column = val;
        else if (arg == "--k") o.cfg.k = integer(val, arg);
        else if (arg == "--min-types") o.cfg.min_types = integer(val, arg);
        else if (arg == "--max-relabels") o.cfg.max_relabels = integer(val, arg);
        else if (arg == "--min-support") o.cfg.min_support = integer(val, arg);
        else if (arg == "--max-motifs") o.cfg.max_motifs = integer(val, arg);
        else if (arg == "--candidate-pool") o.candidate_pool = integer(val, arg);
        else if (arg == "--permutations") o.cfg.permutations = integer(val, arg);
        else if (arg == "--threads") o.cfg.threads = integer(val, arg);
        else if (arg == "--seed") {
            int seed = integer(val, arg);
            if (seed < 0) throw std::runtime_error("seed must be nonnegative");
            o.cfg.seed = static_cast<std::uint64_t>(seed);
        }
        else if (arg == "--max-radius") o.cfg.max_radius = number(val, arg);
        else if (arg == "--max-dominance") o.cfg.max_dominance = number(val, arg);
        else if (arg == "--contact-radius") o.cfg.contact_radius = number(val, arg);
        else if (arg == "--min-cross-edge-fraction") o.cfg.min_cross_edge_fraction = number(val, arg);
        else if (arg == "--region-size") o.cfg.region_size = number(val, arg);
        else if (arg == "--block-size") o.block_size = number(val, arg);
        else if (arg == "--alpha") o.cfg.alpha = number(val, arg);
        else throw std::runtime_error("Unknown option: " + arg);
    }
    const auto& c = o.cfg;
    if (o.input.empty() || o.output.empty()) throw std::runtime_error("--input and --output are required");
    if (c.k < 2 || c.k > 1000 || c.min_types < 2 || c.min_types > c.k)
        throw std::runtime_error("Require 2 <= min-types <= k <= 1000");
    if (c.max_radius < 0 || c.contact_radius < 0 || c.region_size <= 0 || o.block_size <= 0)
        throw std::runtime_error("Radii must be nonnegative and tile sizes positive");
    if (c.max_dominance <= 0 || c.max_dominance > 1 || c.min_cross_edge_fraction < 0 || c.min_cross_edge_fraction > 1)
        throw std::runtime_error("Invalid dominance or cross-edge fraction");
    if (c.min_cross_edge_fraction > 0 && c.contact_radius <= 0)
        throw std::runtime_error("Edge filtering requires --contact-radius > 0");
    if (c.max_relabels < 0 || c.max_relabels > 2 || c.min_support < 1 || c.max_motifs < 1 ||
        c.permutations < 1 || c.threads < 1 || o.candidate_pool < 0 || c.alpha <= 0 || c.alpha >= 1)
        throw std::runtime_error("Invalid tolerance, support, pool, permutations, threads or alpha");
    if (o.candidate_pool && o.candidate_pool < c.max_motifs)
        throw std::runtime_error("candidate-pool must be zero or >= max-motifs");
    if (o.null_model != "global" && o.null_model != "stratified" && o.null_model != "blocks")
        throw std::runtime_error("null must be global, stratified or blocks");
    if (fs::exists(o.output) && !fs::is_directory(o.output))
        throw std::runtime_error("Output must be a directory");
    if (fs::exists(o.output) && !fs::is_empty(o.output) && !o.overwrite)
        throw std::runtime_error("Output directory is not empty; choose a new path or --overwrite");
    return o;
}

// Read RFC4180 quoted fields, including escaped quotes and embedded newlines.
bool csv_row(std::istream& in, std::vector<std::string>& row) {
    row.clear(); std::string field; bool quoted = false, started = false;
    for (;;) {
        const int z = in.get();
        if (z == EOF) {
            if (quoted) throw std::runtime_error("Unterminated quoted CSV field");
            if (!started) return false;
            row.push_back(field); return true;
        }
        started = true; const char c = static_cast<char>(z);
        if (quoted) {
            if (c == '"') {
                if (in.peek() == '"') { in.get(); field += '"'; }
                else quoted = false;
            } else field += c;
        } else if (c == '"' && field.empty()) quoted = true;
        else if (c == ',') { row.push_back(field); field.clear(); }
        else if (c == '\n') { row.push_back(field); return true; }
        else if (c == '\r') {
            if (in.peek() == '\n') in.get();
            row.push_back(field); return true;
        } else field += c;
    }
}
std::map<std::string,int> header(std::istream& in) {
    std::vector<std::string> row;
    if (!csv_row(in,row)) throw std::runtime_error("Empty CSV file");
    if (!row.empty() && row[0].compare(0,3,"\xef\xbb\xbf") == 0) row[0].erase(0,3);
    std::map<std::string,int> out;
    for (int i = 0; i < static_cast<int>(row.size()); ++i)
        if (!out.emplace(trim(row[i]),i).second) throw std::runtime_error("Duplicate CSV column " + row[i]);
    return out;
}
int column(const std::map<std::string,int>& h, const std::string& name, bool required = true) {
    auto it = h.find(name);
    if (it != h.end()) return it->second;
    if (required) throw std::runtime_error("Missing CSV column: " + name);
    return -1;
}
int intern(std::vector<std::string>& names, const std::string& name) {
    auto it = std::find(names.begin(), names.end(), name);
    if (it != names.end()) return static_cast<int>(it - names.begin());
    names.push_back(name); return static_cast<int>(names.size() - 1);
}
Dataset load(const Options& o, int& dropped) {
    std::ifstream f(o.input);
    if (!f) throw std::runtime_error("Cannot read " + o.input);
    const auto h = header(f);
    const int x = column(h,o.x_column), y = column(h,o.y_column), t = column(h,o.type_column);
    const int id = column(h,"Cell_ID",false), s = column(h,"Sample",false);
    const int st = column(h,"Stratum",o.null_model == "stratified"), cp = column(h,"Component",false);
    Dataset d; std::vector<std::string> row; int line = 1;
    std::set<std::pair<int,std::string>> seen;
    while (csv_row(f,row)) {
        ++line;
        if (row.size() == 1 && trim(row[0]).empty()) continue;
        if (row.size() != h.size()) throw std::runtime_error("Wrong CSV field count at record " + std::to_string(line));
        std::string label = trim(row[t]);
        if (lower(label) == "unclassified") { ++dropped; continue; }
        if (label.empty()) throw std::runtime_error("Empty cell type at record " + std::to_string(line));
        Cell c;
        c.x = number(trim(row[x]),"x at record " + std::to_string(line));
        c.y = number(trim(row[y]),"y at record " + std::to_string(line));
        // Keep coordinate magnitudes within reliable integer tile conversion.
        if (std::abs(c.x) > 1e12 || std::abs(c.y) > 1e12)
            throw std::runtime_error("Coordinate magnitude exceeds 1e12");
        c.type = intern(d.types,label);
        c.sample = intern(d.samples,s < 0 ? fs::path(o.input).stem().string() : trim(row[s]));
        c.stratum = intern(d.strata,st < 0 ? "all" : trim(row[st]));
        c.component = intern(d.components,cp < 0 ? "all" : trim(row[cp]));
        c.id = id < 0 ? std::to_string(line-2) : trim(row[id]);
        if (c.id.empty() || d.samples[c.sample].empty() || d.strata[c.stratum].empty() || d.components[c.component].empty())
            throw std::runtime_error("Empty cell ID or grouping value at record " + std::to_string(line));
        if (!seen.insert({c.sample,c.id}).second) throw std::runtime_error("Duplicate cell ID within sample: " + c.id);
        d.cells.push_back(c);
    }
    if (d.cells.empty()) throw std::runtime_error("No classified cells found");
    // Lexical label order makes the count-vector vocabulary independent of first occurrence.
    auto sorted = d.types; std::sort(sorted.begin(),sorted.end());
    for (auto& c : d.cells)
        c.type = static_cast<int>(std::lower_bound(sorted.begin(),sorted.end(),d.types[c.type])-sorted.begin());
    d.types = std::move(sorted);
    return d;
}
struct TemplateSet { std::vector<std::string> ids; std::vector<std::vector<int>> counts; };
TemplateSet load_templates(const Options& o, Dataset& d) {
    TemplateSet out;
    if (o.templates.empty()) return out;
    std::ifstream f(o.templates);
    if (!f) throw std::runtime_error("Cannot read templates " + o.templates);
    auto h = header(f);
    int mi=column(h,"motif_id"), ty=column(h,"cell_type"), co=column(h,"count");
    std::map<std::string,std::map<std::string,int>> rows;
    std::vector<std::string> row;
    while (csv_row(f,row)) {
        if (row.size()==1 && trim(row[0]).empty()) continue;
        if (row.size()!=h.size()) throw std::runtime_error("Wrong template field count");
        const auto id=trim(row[mi]), type=trim(row[ty]);
        int count=integer(trim(row[co]),"template count");
        if (id.empty() || type.empty() || count<0 || count>o.cfg.k) throw std::runtime_error("Invalid template entry");
        if (!rows[id].emplace(type,count).second) throw std::runtime_error("Duplicate template motif/type entry");
        intern(d.types,type); // Absent validation types have zero observed counts.
    }
    for (const auto& [id, entries] : rows) {
        std::vector<int> counts(d.types.size(),0);
        for (const auto& [type,count] : entries) counts[intern(d.types,type)] = count;
        if (std::accumulate(counts.begin(),counts.end(),0) != o.cfg.k || !is_eligible_counts(counts,o.cfg))
            throw std::runtime_error("Template " + id + " violates k or composition eligibility");
        out.ids.push_back(id); out.counts.push_back(std::move(counts));
    }
    if (out.ids.empty()) throw std::runtime_error("Template file contains no motifs");
    return out;
}
std::vector<std::vector<int>> shuffle_groups(const Dataset& d, const Options& o) {
    using Key=std::tuple<int,int,int,long long,long long>;
    std::map<Key,std::vector<int>> groups;
    for (int i=0;i<static_cast<int>(d.cells.size());++i) {
        const auto& c=d.cells[i];
        int st=o.null_model=="global" ? 0 : c.stratum;
        long long bx=0,by=0;
        if (o.null_model=="blocks") {
            const double fx=std::floor(c.x/o.block_size),fy=std::floor(c.y/o.block_size);
            if (std::abs(fx)>9e18 || std::abs(fy)>9e18) throw std::runtime_error("block-size is too small for coordinates");
            bx=static_cast<long long>(fx); by=static_cast<long long>(fy);
        }
        groups[{c.sample,c.component,st,bx,by}].push_back(i);
    }
    std::vector<std::vector<int>> out;
    for (auto& [key,idx] : groups) out.push_back(std::move(idx));
    return out;
}
std::uint64_t splitmix(std::uint64_t x) {
    x+=0x9e3779b97f4a7c15ULL; x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
    x=(x^(x>>27))*0x94d049bb133111ebULL; return x^(x>>31);
}
std::string quote(const std::string& s) {
    std::string out="\"";
    for (char c:s) { if(c=='"') out+='"'; out+=c; }
    return out+'"';
}
std::string json(const std::string& s) {
    std::ostringstream out; out << '"';
    for (unsigned char c:s) {
        if(c=='"'||c=='\\') out << '\\' << c;
        else if(c=='\n') out << "\\n";
        else if(c=='\r') out << "\\r";
        else if(c=='\t') out << "\\t";
        else if(c<32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c) << std::dec;
        else out << c;
    }
    out << '"'; return out.str();
}
std::ofstream output(const Options& o,const std::string& name) {
    std::ofstream f(fs::path(o.output)/name);
    if(!f) throw std::runtime_error("Cannot write " + name);
    f << std::setprecision(12); return f;
}
std::string composition(const Dataset& d,const std::vector<int>& counts) {
    std::ostringstream out; bool first=true;
    for(std::size_t t=0;t<counts.size();++t) if(counts[t]) {
        if(!first) out << "; ";
        first=false;
        out << d.types[t] << ':' << counts[t];
    }
    return out.str();
}
bool same_type_set(const std::vector<int>& a,const std::vector<int>& b) {
    for(std::size_t t=0;t<a.size();++t) if((a[t]>0)!=(b[t]>0)) return false;
    return true;
}
double pvalue(int observed,const std::vector<int>& null) {
    return (1.0+std::count_if(null.begin(),null.end(),[&](int n){return n>=observed;}))/(null.size()+1.0);
}
int eligible_count(const std::vector<NeighborhoodStats>& stats) {
    return static_cast<int>(std::count_if(stats.begin(),stats.end(),[](const auto& s){return s.eligible;}));
}
} // namespace

int main(int argc,char** argv) {
    try {
        const auto started=std::chrono::steady_clock::now();
        const auto o=parse(argc,argv); const auto& c=o.cfg;
        int dropped=0; auto d=load(o,dropped); const bool validation=!o.templates.empty();
        auto fixed=load_templates(o,d);
        std::vector<int> labels; for(const auto& cell:d.cells) labels.push_back(cell.type);
        auto groups=shuffle_groups(d,o);
        int exchangeable_cells=0,mixed_strata=0;
        for(const auto& g:groups) {
            bool mixed=false; for(int i:g) if(labels[i]!=labels[g[0]]) { mixed=true; break; }
            if(mixed) { exchangeable_cells+=static_cast<int>(g.size()); ++mixed_strata; }
        }
        std::cerr << "Building " << c.k << "-cell neighborhoods for " << d.cells.size() << " cells\n";
        const auto geom=build_geometry(d,c);
        const auto observed_stats=characterize(d,geom,labels,c);
        Config search_config=c;
        search_config.min_support=1; // Preserve the full null distribution below reporting floor.
        auto observed=validation ? score_templates(d,geom,observed_stats,fixed.counts,c)
                                 : mine(d,geom,observed_stats,search_config);
        const int searched_count=static_cast<int>(observed.size());
        const int observed_max=observed.empty()?0:std::max_element(observed.begin(),observed.end(),
            [](const auto& a,const auto& b){return a.support<b.support;})->support;
        if(!validation) observed.erase(std::remove_if(observed.begin(),observed.end(),
            [&](const auto& motif){return motif.support<c.min_support;}),observed.end());
        const int candidate_count=static_cast<int>(observed.size());
        if(!validation && o.candidate_pool && static_cast<int>(observed.size())>o.candidate_pool)
            observed.resize(o.candidate_pool);
        std::vector<std::vector<int>> templates;
        for(const auto& motif:observed) templates.push_back(motif.counts);
        const int m=static_cast<int>(templates.size()),B=c.permutations;
        std::vector<int> maxima(B),null_eligible(B);
        std::vector<std::vector<int>> motif_null(m,std::vector<int>(B));
        std::atomic<int> next{0},done{0}; std::mutex lock; std::exception_ptr failure;
        std::atomic<bool> failed{false};
        std::cerr << geom.neighborhoods.size() << " distinct neighborhoods; " << eligible_count(observed_stats)
                  << " mixed eligible; " << candidate_count << " recurring candidates; " << m << " templates in enrichment pool\n";
        const int workers=std::min(c.threads,B);
        std::vector<std::thread> threads;
        for(int worker=0;worker<workers;++worker) threads.emplace_back([&] {
            try {
                while(!failed) {
                    const int b=next.fetch_add(1); if(b>=B) break;
                    std::mt19937_64 rng(splitmix(c.seed+static_cast<std::uint64_t>(b)));
                    auto shuffled=labels;
                    for(const auto& g:groups) {
                        // Fisher-Yates within each fixed sample/stratum, including all cells.
                        for(std::size_t j=g.size();j>1;--j) {
                            std::uniform_int_distribution<std::size_t> draw(0,j-1);
                            std::swap(shuffled[g[j-1]],shuffled[g[draw(rng)]]);
                        }
                    }
                    const auto stats=characterize(d,geom,shuffled,c);
                    null_eligible[b]=eligible_count(stats);
                    const auto scored=score_templates(d,geom,stats,templates,c);
                    for(int j=0;j<m;++j) motif_null[j][b]=scored[j].support;
                    if(validation) {
                        for(const auto& motif:scored) maxima[b]=std::max(maxima[b],motif.support);
                    } else {
                        const auto mined=mine(d,geom,stats,search_config);
                        if(!mined.empty()) maxima[b]=mined.front().support;
                    }
                    const int count=done.fetch_add(1)+1;
                    if(count%25==0 || count==B) {
                        std::lock_guard<std::mutex> guard(lock);
                        std::cerr << "Permutation " << count << '/' << B << '\n';
                    }
                }
            } catch(...) {
                std::lock_guard<std::mutex> guard(lock); if(!failure) failure=std::current_exception(); failed=true;
            }
        });
        for(auto& t:threads) t.join();
        if(failure) std::rethrow_exception(failure);
        const double global_p=pvalue(observed_max,maxima);
        std::vector<double> means(m),p(m),adjusted(m),excess(m);
        for(int j=0;j<m;++j) {
            means[j]=std::accumulate(motif_null[j].begin(),motif_null[j].end(),0.0)/B;
            excess[j]=observed[j].support-means[j];
            p[j]=pvalue(observed[j].support,motif_null[j]);
        }
        std::vector<int> order(m); std::iota(order.begin(),order.end(),0);
        if(validation) {
            std::sort(order.begin(),order.end(),[&](int a,int b){return p[a]<p[b] || (p[a]==p[b] && a<b);});
            double running=0;
            for(int r=0;r<m;++r) { running=std::max(running,std::min(1.0,(m-r)*p[order[r]])); adjusted[order[r]]=running; }
            std::iota(order.begin(),order.end(),0);
        } else {
            std::sort(order.begin(),order.end(),[&](int a,int b) {
                if(excess[a]!=excess[b]) return excess[a]>excess[b];
                if(observed[a].support!=observed[b].support) return observed[a].support>observed[b].support;
                return observed[a].counts<observed[b].counts;
            });
        }
        std::vector<int> selected;
        for(int j:order) {
            bool redundant=false;
            if(!validation) for(int z:selected)
                if((!c.same_types || same_type_set(templates[j],templates[z])) &&
                    relabel_distance(templates[j],templates[z])<=2*c.max_relabels) {redundant=true;break;}
            if(!redundant) selected.push_back(j);
            if(!validation && static_cast<int>(selected.size())==c.max_motifs) break;
        }
        std::vector<std::string> ids(m);
        for(std::size_t r=0;r<selected.size();++r) {
            auto j=selected[r]; std::ostringstream name; name << 'M' << std::setw(3) << std::setfill('0') << r+1;
            ids[j]=validation?fixed.ids[j]:name.str();
        }
        fs::create_directories(o.output);
        {
            auto f=output(o,"cells.csv"); f << "cell_index,cell_id,x,y,cell_type,sample,stratum,component\n";
            for(std::size_t i=0;i<d.cells.size();++i) {const auto& cell=d.cells[i];
                f << i << ',' << quote(cell.id) << ',' << cell.x << ',' << cell.y << ',' << quote(d.types[cell.type]) << ','
                  << quote(d.samples[cell.sample]) << ',' << quote(d.strata[cell.stratum]) << ',' << quote(d.components[cell.component]) << '\n';}
        }
        {
            auto f=output(o,"neighborhoods.csv");
            f << "neighborhood_id,center_index,x,y,sample,size,radius,n_types,dominance,entropy,effective_types,cross_edge_fraction,eligible\n";
            for(std::size_t i=0;i<geom.neighborhoods.size();++i) {
                const auto& n=geom.neighborhoods[i];const auto& s=observed_stats[i];const auto& cell=d.cells[n.center];
                f << i << ',' << n.center << ',' << cell.x << ',' << cell.y << ',' << quote(d.samples[cell.sample]) << ',' << n.members.size()
                  << ',' << n.radius << ',' << s.n_types << ',' << s.dominance << ',' << s.entropy << ',' << s.effective_types
                  << ',' << s.cross_edge_fraction << ',' << (s.eligible?"true":"false") << '\n';
            }
        }
        const std::string scope=validation ? "fixed_template_random_label_test_requires_independent_specification"
            : "complete_random_label_null_search_calibration_not_motifwise_FDR";
        {
            auto f=output(o,"motifs.csv"), comp=output(o,"motif_composition.csv"),occ=output(o,"occurrences.csv"),mem=output(o,"members.csv");
            f << "motif_id,composition,support,raw_matches,n_regions,n_samples,expected_support,fold_enrichment,excess_support,p_value,p_value_holm,p_search,significant,inference_scope\n";
            comp << "motif_id,cell_type,count\n";
            occ << "motif_id,occurrence_id,neighborhood_id,center_index,x,y,sample,radius\n";
            mem << "motif_id,occurrence_id,cell_index\n";
            for(int j:selected) {
                const auto& motif=observed[j];const double search_p=pvalue(motif.support,maxima);
                const bool significant=(validation?adjusted[j]<=c.alpha:search_p<=c.alpha) && motif.support>=c.min_support;
                f << quote(ids[j]) << ',' << quote(composition(d,motif.counts)) << ',' << motif.support << ',' << motif.raw_matches << ','
                  << motif.regions << ',' << motif.samples << ',' << means[j] << ',';
                // Finite 0.5 pseudocount ratio, recorded explicitly in run metadata.
                f << (motif.support+0.5)/(means[j]+0.5) << ',' << excess[j] << ',';
                if(validation) f << p[j] << ',' << adjusted[j] << ',';
                else f << ",,";
                if(!validation) f << search_p;
                f << ',' << (significant?"true":"false") << ',' << quote(scope) << '\n';
                for(std::size_t t=0;t<d.types.size();++t) comp << quote(ids[j]) << ',' << quote(d.types[t]) << ',' << motif.counts[t] << '\n';
                for(std::size_t r=0;r<motif.occurrences.size();++r) {
                    int ni=motif.occurrences[r];const auto& n=geom.neighborhoods[ni];const auto& cell=d.cells[n.center];
                    occ << quote(ids[j]) << ',' << r+1 << ',' << ni << ',' << n.center << ',' << cell.x << ',' << cell.y << ',' << quote(d.samples[cell.sample]) << ',' << n.radius << '\n';
                    for(int member:n.members) mem << quote(ids[j]) << ',' << r+1 << ',' << member << '\n';
                }
            }
        }
        {
            auto f=output(o,"null.csv"), mn=output(o,"motif_null.csv");
            f << "permutation,max_support,eligible_neighborhoods\n";
            mn << "permutation,motif_id,support\n";
            for(int b=0;b<B;++b) {
                f << b+1 << ',' << maxima[b] << ',' << null_eligible[b] << '\n';
                for(int j:selected) mn << b+1 << ',' << quote(ids[j]) << ',' << motif_null[j][b] << '\n';
            }
        }
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        {
            auto f=output(o,"run.json");
            f << "{\n  \"schema_version\": 1,\n  \"version\": \"1.0.0\",\n  \"mode\": " << json(validation?"validate":"discovery")
              << ",\n  \"input\": " << json(fs::absolute(o.input).string()) << ",\n  \"templates\": " << json(o.templates)
              << ",\n  \"null_model\": " << json(o.null_model)
              << ",\n  \"null_assumption\": " << json("Labels exchangeable within fixed sample/component/strata; within-stratum homotypic clustering is not preserved.")
              << ",\n  \"inference_scope\": " << json(scope)
              << ",\n  \"support_definition\": \"compact-first greedy cell-disjoint occurrences; not independent biological replicates\""
              << ",\n  \"rank_by\": " << json(validation?"template_id":"excess_support_in_support_screened_pool")
              << ",\n  \"enrichment_interpretation\": " << json(validation?"fixed-template null expectation":"descriptive; observed-template selection and ranking reused data")
              << ",\n  \"fold_pseudocount\": 0.5,\n  \"cells\": " << d.cells.size() << ",\n  \"dropped_unclassified\": " << dropped
              << ",\n  \"n_neighborhoods\": " << geom.neighborhoods.size() << ",\n  \"n_eligible\": " << eligible_count(observed_stats)
              << ",\n  \"rejected_radius\": " << geom.rejected_radius << ",\n  \"duplicate_sets\": " << geom.duplicate_sets
              << ",\n  \"n_searched_templates\": " << searched_count
              << ",\n  \"n_candidates\": " << candidate_count << ",\n  \"candidate_pool\": " << o.candidate_pool
              << ",\n  \"n_scored_templates\": " << m << ",\n  \"n_reported_motifs\": " << selected.size()
              << ",\n  \"observed_max_support\": " << observed_max << ",\n  \"global_p_value\": ";
            if(validation) f << "null"; else f << global_p;
            f << ",\n  \"k\": " << c.k << ",\n  \"max_radius\": " << c.max_radius << ",\n  \"min_types\": " << c.min_types
              << ",\n  \"max_dominance\": " << c.max_dominance << ",\n  \"max_relabels\": " << c.max_relabels
              << ",\n  \"same_types\": " << (c.same_types?"true":"false")
              << ",\n  \"min_support\": " << c.min_support << ",\n  \"max_motifs\": " << c.max_motifs
              << ",\n  \"contact_radius\": " << c.contact_radius << ",\n  \"min_cross_edge_fraction\": " << c.min_cross_edge_fraction
              << ",\n  \"region_size\": " << c.region_size << ",\n  \"block_size\": " << o.block_size
              << ",\n  \"n_shuffle_strata\": " << groups.size() << ",\n  \"n_mixed_shuffle_strata\": " << mixed_strata
              << ",\n  \"exchangeable_cells\": " << exchangeable_cells << ",\n  \"permutations\": " << B
              << ",\n  \"minimum_p_value\": " << 1.0/(B+1) << ",\n  \"alpha\": " << c.alpha
              << ",\n  \"seed\": " << c.seed << ",\n  \"threads\": " << workers << ",\n  \"elapsed_seconds\": " << seconds << "\n}\n";
        }
        std::cout << "Results: " << fs::absolute(o.output) << "\nReported motifs: " << selected.size() << '\n';
        if(!validation) std::cout << "Global random-label search p-value: " << global_p << '\n';
        std::cout << "Elapsed seconds: " << seconds << '\n';
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n'; return 1;
    }
}
