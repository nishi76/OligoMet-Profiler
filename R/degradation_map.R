# =============================================================================
# degradation_map.R
# Position-resolved degradation map: every identified metabolite is placed on
# the parent sequence, colored by degradation class, and annotated by
# confidence, across the study's timepoints and matrices. Modeled on the
# Sprinzl-position modification map of Shen et al. (OligoFragger, 2026, Fig 3d),
# adapted from endogenous-RNA modification discovery to therapeutic-oligo
# metID/DMPK: rows are timepoints, columns are residue (and linkage) positions,
# cell color is the degradation class, cell opacity scales with abundance
# relative to the full-length parent.
#
# degradation_position_model() is pure base R and is the testable core: it
# turns the batch match table + metabolite library + sample metadata into a
# tidy set of cells. degradation_map_html() renders those cells as a
# self-contained HTML/JS widget for the Shiny app and the study report.
# =============================================================================

# Canonical cell-class keys consumed by the widget (color + glyph live in JS).
# Base-modification grid: the three position-resolved in-place mods.
.DM_MOD_CLASS <- c(depurination = "aB", oxidation = "ox", deamination = "dA")
# Backbone cleavage lane, keyed by metabolite kind family.
.DM_CLEAVE_CLASS <- c(exo_3p = "e3", exo_5p = "e5",
                      endo_5frag = "en", endo_3frag = "en",
                      endo_5frag_p = "en", endo_3frag_p = "en")

# Residue span [start, end] (1-indexed, on the parent) a metabolite covers, so
# a position can be marked "covered but unmodified" when a detected species
# spans it. Derived from kind + site + parent length; NULL when not placeable.
.dm_met_span <- function(kinds, site, parent_n) {
  k <- kinds[1]
  s <- suppressWarnings(as.integer(site))
  if (k %in% c("parent", "conj_catab")) return(c(1L, parent_n))
  if (is.na(s)) return(NULL)
  if (k == "exo_3p") return(c(1L, max(1L, parent_n - s)))
  if (k == "exo_5p") return(c(min(parent_n, s + 1L), parent_n))
  if (k %in% c("endo_5frag", "endo_5frag_p")) return(c(1L, s))
  if (k %in% c("endo_3frag", "endo_3frag_p")) return(c(min(parent_n, s + 1L), parent_n))
  NULL
}

# Parse a timepoint label to hours for ordering. Accepts "predose"/"pre"/"0",
# "4h"/"4 h", "30 d"/"30 days", "1 week", bare numbers (assumed hours). NA when
# unparseable, which sorts last in appearance order.
.dm_parse_hours <- function(x) {
  v <- tolower(trimws(as.character(x)))
  if (is.na(v) || !nzchar(v)) return(NA_real_)
  if (v %in% c("predose", "pre", "pre-dose", "baseline", "t0")) return(0)
  num <- suppressWarnings(as.numeric(sub("^\\s*([0-9.]+).*$", "\\1", v)))
  if (is.na(num)) return(NA_real_)
  if (grepl("day|\\bd\\b|\\bd$", v)) return(num * 24)
  if (grepl("week|\\bwk\\b|\\bw$", v)) return(num * 24 * 7)
  if (grepl("min|\\bm$", v) && !grepl("month", v)) return(num / 60)
  num  # hours, or a bare number
}

.dm_order_levels <- function(labels) {
  u <- unique(labels[!is.na(labels) & nzchar(labels)])
  if (length(u) == 0) return(character(0))
  hrs <- vapply(u, .dm_parse_hours, numeric(1))
  u[order(ifelse(is.na(hrs), Inf, hrs), seq_along(u))]
}

.dm_signal_col <- function(matches, signal_col) {
  if (signal_col %in% names(matches) && any(!is.na(matches[[signal_col]]))) return(signal_col)
  if ("area" %in% names(matches) && any(!is.na(matches$area))) return("area")
  if ("intensity" %in% names(matches) && any(!is.na(matches$intensity))) return("intensity")
  NA_character_
}

#' Build the position-resolved degradation-map model
#'
#' Turns the batch match table into a tidy set of cells for the degradation
#' map: in-place base modifications placed on their residue, backbone
#' cleavage placed on its linkage, and the GalNAc conjugate ladder, each
#' with abundance relative to the full-length parent, per timepoint and
#' matrix. Every piece of information comes from objects the pipeline
#' already produces (`annotate_metabolites_batch()` matches, the metabolite
#' library, and the sample metadata); this function only reshapes them.
#'
#' @param matches A batch match data.frame (`ms1_matches` from
#'   [annotate_metabolites_batch()]): needs `met_id`, `kind`, `sample`, a
#'   signal column, and optionally `mod`, `mod_site`, `ppm_error`,
#'   `ms2_score`, `n`.
#' @param mets The metabolite library (see [generate_metabolites()]), used
#'   for each species' residue span and the parent length.
#' @param spec The parent `oligo_spec` (sequence, sugars for the DNA gap).
#' @param sample_meta Optional data.frame keyed by `sample` with
#'   `timepoint` and/or `matrix` columns (see [read_sample_info()]). When
#'   absent, all samples fall under one timepoint/matrix.
#' @param oligo_name Parent display name.
#' @param signal_col Preferred signal column (`"area"`, falls back to
#'   `"intensity"`).
#' @param min_pct Drop cells below this percent of parent (default 0.1).
#' @param ms2_cutoff `ms2_score` at or above which a cell is marked
#'   MS2-confirmed (default 0).
#' @return A list with `oligo`, `seq`, `gap`, `times`, `matrices`, and the
#'   cell frames `base`, `cleave`, `conj`, `covered`, `parent` (parent
#'   percent remaining vs the earliest timepoint). Empty frames when there
#'   is nothing to place.
#' @seealso [degradation_map_html()], [degradation_by_condition()]
#' @export
degradation_position_model <- function(matches, mets, spec, sample_meta = NULL,
                                       oligo_name = "parent", signal_col = "area",
                                       min_pct = 0.1, ms2_cutoff = 0) {
  parent_n <- spec$n %||% length(spec$bases)
  gap <- NULL
  if (!is.null(spec$sugars)) {
    d <- which(spec$sugars == "d")
    if (length(d) > 0 && !identical(spec$sugars[1], "d") &&
        !identical(spec$sugars[parent_n], "d")) gap <- c(min(d), max(d))
  }
  empty <- function() data.frame()
  model <- list(oligo = oligo_name, seq = as.character(spec$bases),
                gap = gap, times = character(0), matrices = character(0),
                base = empty(), cleave = empty(), conj = empty(),
                covered = empty(), parent = empty())
  if (is.null(matches) || nrow(matches) == 0) return(model)

  sc_use <- .dm_signal_col(matches, signal_col)
  if (is.na(sc_use)) return(model)
  m <- matches
  m$.sig <- suppressWarnings(as.numeric(m[[sc_use]]))
  m <- m[!is.na(m$.sig) & m$.sig > 0, , drop = FALSE]
  if (nrow(m) == 0) return(model)

  # kind lookup (kind_all carries every route after dedupe)
  kind_by_id <- setNames(lapply(mets, function(x) x$kind_all %||% x$kind),
                         vapply(mets, function(x) x$id, ""))
  site_by_id <- setNames(lapply(mets, function(x) x$site), vapply(mets, function(x) x$id, ""))

  # sample -> timepoint / matrix
  tp <- rep("all", nrow(m)); mx <- rep("all", nrow(m))
  if (!is.null(sample_meta) && "sample" %in% names(sample_meta)) {
    idx <- match(m$sample, sample_meta$sample)
    if ("timepoint" %in% names(sample_meta))
      tp <- ifelse(is.na(idx), "all", as.character(sample_meta$timepoint[idx]))
    if ("matrix" %in% names(sample_meta))
      mx <- ifelse(is.na(idx), "all", as.character(sample_meta$matrix[idx]))
    tp[is.na(tp) | !nzchar(tp)] <- "all"; mx[is.na(mx) | !nzchar(mx)] <- "all"
  }
  m$.tp <- tp; m$.mx <- mx

  # per-sample parent signal (max over parent rows in that sample)
  is_parent <- vapply(m$kind, function(k) identical(k, "parent"), logical(1))
  parent_sig <- tapply(m$.sig[is_parent], m$sample[is_parent], max)
  m$.psig <- as.numeric(parent_sig[m$sample])
  m$.pct <- ifelse(is.na(m$.psig) | m$.psig <= 0, NA_real_, 100 * m$.sig / m$.psig)

  mod <- if ("mod" %in% names(m)) as.character(m$mod) else rep("none", nrow(m))
  mod[is.na(mod)] <- "none"
  msite <- if ("mod_site" %in% names(m)) suppressWarnings(as.integer(m$mod_site)) else rep(NA_integer_, nrow(m))
  ppm <- if ("ppm_error" %in% names(m)) suppressWarnings(as.numeric(m$ppm_error)) else rep(NA_real_, nrow(m))
  ms2 <- if ("ms2_score" %in% names(m)) suppressWarnings(as.numeric(m$ms2_score)) else rep(NA_real_, nrow(m))
  ms2ok <- !is.na(ms2) & ms2 >= ms2_cutoff

  base_rows <- list(); cleave_rows <- list(); conj_rows <- list(); cov_rows <- list()
  add <- function(lst, row) { lst[[length(lst) + 1L]] <- row; lst }

  for (i in seq_len(nrow(m))) {
    pct <- m$.pct[i]; if (is.na(pct)) next
    kinds <- kind_by_id[[m$met_id[i]]] %||% m$kind[i]
    tcol <- m$.tp[i]; xcol <- m$.mx[i]
    # coverage: every position this species spans, when detected
    sp <- .dm_met_span(kinds, site_by_id[[m$met_id[i]]] %||% m$n[i], parent_n)
    if (!is.null(sp)) for (p in sp[1]:sp[2])
      cov_rows <- add(cov_rows, data.frame(matrix = xcol, time = tcol, pos = p, stringsAsFactors = FALSE))
    # 1) in-place base modification -> grid cell at mod_site
    if (mod[i] %in% names(.DM_MOD_CLASS) && !is.na(msite[i])) {
      base_rows <- add(base_rows, data.frame(
        matrix = xcol, time = tcol, pos = msite[i], cls = unname(.DM_MOD_CLASS[mod[i]]),
        pct = pct, ms2 = ms2ok[i], cyc = FALSE, ppm = ppm[i],
        species = m$met_name[i], stringsAsFactors = FALSE))
      next
    }
    # 2) cyclic phosphate -> cleavage lane flag at the species' 3' linkage
    cyc <- identical(mod[i], "cyclic_phosphate")
    # 3) backbone cleavage -> lane cell at the cleavage linkage (metabolite site)
    kf <- intersect(kinds, names(.DM_CLEAVE_CLASS))
    if (length(kf) > 0) {
      link <- suppressWarnings(as.integer(site_by_id[[m$met_id[i]]]))
      if (is.na(link)) link <- suppressWarnings(as.integer(m$n[i]))
      cleave_rows <- add(cleave_rows, data.frame(
        matrix = xcol, time = tcol, link = link, cls = unname(.DM_CLEAVE_CLASS[kf[1]]),
        pct = pct, ms2 = ms2ok[i], cyc = cyc, species = m$met_name[i], stringsAsFactors = FALSE))
      next
    }
    # 4) conjugate catabolism -> conjugate track
    if ("conj_catab" %in% kinds) {
      conj_rows <- add(conj_rows, data.frame(
        matrix = xcol, time = tcol, cls = "gn", pct = pct, ms2 = ms2ok[i],
        label = m$met_name[i], species = m$met_name[i], stringsAsFactors = FALSE))
    }
  }

  bind <- function(lst) if (length(lst) == 0) data.frame() else do.call(rbind, lst)
  base <- bind(base_rows); cleave <- bind(cleave_rows)
  conj <- bind(conj_rows); cov <- bind(cov_rows)

  # aggregate to one cell per (matrix, time, position/linkage, class): max %
  agg_cell <- function(df, poscol) {
    if (nrow(df) == 0) return(df)
    key <- paste(df$matrix, df$time, df[[poscol]], df$cls, sep = "\r")
    ord <- order(key, -df$pct)
    df <- df[ord, , drop = FALSE]
    df[!duplicated(paste(df$matrix, df$time, df[[poscol]], df$cls, sep = "\r")), , drop = FALSE]
  }
  base <- agg_cell(base, "pos"); cleave <- agg_cell(cleave, "link")
  if (nrow(base) > 0) base <- base[base$pct >= min_pct, , drop = FALSE]
  if (nrow(cleave) > 0) cleave <- cleave[cleave$pct >= min_pct, , drop = FALSE]
  if (nrow(conj) > 0) {
    key <- paste(conj$matrix, conj$time, conj$label, sep = "\r")
    ord <- order(key, -conj$pct); conj <- conj[ord, , drop = FALSE]
    conj <- conj[!duplicated(paste(conj$matrix, conj$time, conj$label, sep = "\r")), , drop = FALSE]
    conj <- conj[conj$pct >= min_pct, , drop = FALSE]
  }
  if (nrow(cov) > 0) cov <- unique(cov)

  # parent percent remaining vs the earliest timepoint, per matrix
  par_rows <- list()
  if (any(is_parent)) {
    pm <- m[is_parent, , drop = FALSE]
    agg <- stats::aggregate(.sig ~ .tp + .mx, data = pm, FUN = max)
    for (xcol in unique(agg$.mx)) {
      sub <- agg[agg$.mx == xcol, , drop = FALSE]
      lv <- .dm_order_levels(sub$.tp)
      base_sig <- sub$.sig[match(lv[1], sub$.tp)]
      if (is.na(base_sig) || base_sig <= 0) next
      for (tt in lv) par_rows[[length(par_rows) + 1L]] <- data.frame(
        matrix = xcol, time = tt, pct = 100 * sub$.sig[match(tt, sub$.tp)] / base_sig,
        stringsAsFactors = FALSE)
    }
  }
  parent <- bind(par_rows)

  all_times <- .dm_order_levels(c(if (nrow(base)) base$time, if (nrow(cleave)) cleave$time,
                                  if (nrow(conj)) conj$time, if (nrow(parent)) parent$time, m$.tp))
  all_mx <- unique(c(if (nrow(base)) base$matrix, if (nrow(cleave)) cleave$matrix,
                     if (nrow(conj)) conj$matrix, m$.mx))
  model$times <- all_times
  model$matrices <- all_mx
  model$base <- base; model$cleave <- cleave
  model$conj <- conj; model$covered <- cov; model$parent <- parent
  model
}

# ---- HTML widget -----------------------------------------------------------
# The widget template: CSS + markup + JS, all scoped under .dmap so it can be
# embedded inside the Shiny app (or the study report) without leaking styles.
# The JS reads a single injected DATA object (see degradation_map_html()).
.DM_WIDGET_CSS <- '
<style>
.dmap{--s:#fcfcfb;--pl:#f6f6f3;--pa:#fff;--ink:#0b0b0b;--ink2:#52514e;--mut:#898781;
 --grid:#e1e0d9;--ring:rgba(11,11,11,.10);
 --ox:#2a78d6;--aB:#eb6834;--dA:#1baf7a;--e3:#6250d6;--e5:#e87ba4;--en:#008300;--gn:#eda100;
 --cov:#d8d7d0;--wing:#b9c6d6;--gap:#e7cdb6;--li:#fff;
 font-family:system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--ink);font-size:13px}
@media (prefers-color-scheme:dark){.dmap{--s:#1a1a19;--pl:#141413;--pa:#1f1f1d;--ink:#f0efec;
 --ink2:#c3c2b7;--mut:#898781;--grid:#2c2c2a;--ring:rgba(255,255,255,.10);
 --ox:#3987e5;--aB:#d95926;--dA:#199e70;--e3:#9085e9;--e5:#d55181;--en:#2fae2f;--gn:#c98500;
 --cov:#38382f;--wing:#3a485c;--gap:#4a3a2c}}
.dmap *{box-sizing:border-box}
.dmap .tabs{display:flex;gap:4px;flex-wrap:wrap;margin:0 0 12px}
.dmap .tab{appearance:none;border:1px solid var(--ring);background:var(--pa);color:var(--ink2);
 font:inherit;font-size:12px;padding:5px 12px;border-radius:999px;cursor:pointer}
.dmap .tab[aria-selected=true]{background:var(--ink);color:var(--s);border-color:var(--ink);font-weight:600}
.dmap .panel{background:var(--s);border:1px solid var(--ring);border-radius:10px;padding:13px;margin-bottom:13px}
.dmap .panel h3{font-size:12px;letter-spacing:.02em;text-transform:uppercase;color:var(--mut);margin:0 0 10px;font-weight:600}
.dmap .scrollx{overflow-x:auto}
.dmap table{border-collapse:separate;border-spacing:0;font-variant-numeric:tabular-nums}
.dmap th,.dmap td{padding:0;text-align:center;vertical-align:middle}
.dmap .rl{position:sticky;left:0;z-index:2;background:var(--s);text-align:right;padding-right:9px;
 font-size:11px;color:var(--ink2);white-space:nowrap;min-width:54px}
.dmap .rl .p{display:block;font-size:10px;color:var(--mut)}
.dmap .ph{font-size:9.5px;color:var(--mut);height:15px;width:26px}
.dmap .sc{width:26px;height:28px;font-size:12px;font-weight:600;position:relative}
.dmap .sc::after{content:"";position:absolute;left:2px;right:2px;bottom:2px;height:3px;border-radius:2px;background:var(--wing)}
.dmap .sc.g::after{background:var(--gap)}
.dmap .c{width:26px;height:26px;position:relative}
.dmap .c .bx{position:absolute;inset:2px;border-radius:5px;display:flex;align-items:center;justify-content:center;
 font-size:9.5px;font-weight:700;color:var(--li);letter-spacing:.02em}
.dmap .c.cov .bx{background:var(--cov);color:var(--mut)}
.dmap .c .d{position:absolute;right:2px;bottom:2px;width:5px;height:5px;border-radius:50%;background:var(--s);box-shadow:0 0 0 1px var(--ring)}
.dmap .c .cp{position:absolute;right:2px;top:0;font-size:8px;color:var(--li)}
.dmap .trk{background:var(--grid);border-radius:5px;height:8px;width:48px;display:inline-block;overflow:hidden;vertical-align:middle}
.dmap .trk>span{display:block;height:8px;background:var(--ox)}
.dmap .legend{display:flex;flex-wrap:wrap;gap:12px 18px}
.dmap .lg{display:flex;flex-direction:column;gap:6px;min-width:140px}
.dmap .lg .t{font-size:10.5px;text-transform:uppercase;letter-spacing:.03em;color:var(--mut);font-weight:600}
.dmap .lr{display:flex;align-items:center;gap:7px;font-size:12px;color:var(--ink2)}
.dmap .sw{width:16px;height:16px;border-radius:4px;display:inline-flex;align-items:center;justify-content:center;
 font-size:8.5px;font-weight:700;color:var(--li);flex:none}
.dmap .sw.cov{background:var(--cov);color:var(--mut)}
.dmap .note{font-size:11.5px;color:var(--mut);margin-top:9px;line-height:1.5}
.dmap .mkdot{width:11px;height:11px;border-radius:50%;background:var(--cov);position:relative;display:inline-block}
.dmap .mkdot::after{content:"";position:absolute;right:0;bottom:0;width:5px;height:5px;border-radius:50%;background:var(--ox);box-shadow:0 0 0 1px var(--s)}
.dmap .tip{position:fixed;z-index:9999;pointer-events:none;opacity:0;transition:opacity .08s;background:var(--pa);
 color:var(--ink);border:1px solid var(--ring);border-radius:8px;padding:7px 9px;font-size:12px;max-width:250px;box-shadow:0 6px 24px rgba(0,0,0,.18)}
.dmap .tip .k{color:var(--mut)}
.dmap .tip hr{border:none;border-top:1px solid var(--grid);margin:5px 0}
.dmap .empty{color:var(--mut);font-size:12.5px;padding:18px 4px}
</style>'

.DM_WIDGET_JS <- '
<script>
(function(){
 var D=__DATA__, root=document.getElementById("__ID__"); if(!root) return;
 var CLS={ox:["--ox","ox","Oxidation"],aB:["--aB","aB","Depurination/abasic"],dA:["--dA","dA","Deamination"],
   e3:["--e3","3\\u2032","3\\u2032-exonuclease"],e5:["--e5","5\\u2032","5\\u2032-exonuclease"],
   en:["--en","en","Endonuclease"],gn:["--gn","Gn","GalNAc loss"]};
 function cvar(v){return getComputedStyle(root).getPropertyValue(v).trim();}
 function rgba(hex,op){var n=hex.replace("#","");return "rgba("+parseInt(n.slice(0,2),16)+","+parseInt(n.slice(2,4),16)+","+parseInt(n.slice(4,6),16)+","+op.toFixed(3)+")";}
 function alpha(p,mx){if(p<=0.05)return 0;return 0.22+0.78*Math.sqrt(Math.min(p,mx)/mx);}
 var tip=root.querySelector(".tip");
 function show(h,x,y){tip.innerHTML=h;tip.style.opacity="1";var w=tip.offsetWidth,hh=tip.offsetHeight,nx=x+14,ny=y+14;
   if(nx+w>innerWidth-8)nx=x-w-14; if(ny+hh>innerHeight-8)ny=y-hh-14; tip.style.left=nx+"px";tip.style.top=ny+"px";}
 function hide(){tip.style.opacity="0";}
 function bind(el,h){el.addEventListener("mousemove",function(e){show(h,e.clientX,e.clientY);});
   el.addEventListener("mouseleave",hide);el.tabIndex=0;
   el.addEventListener("focus",function(){var r=el.getBoundingClientRect();show(h,r.right,r.bottom);});el.addEventListener("blur",hide);}
 var N=D.seq.length;
 function inGap(p){return D.gap&&p>=D.gap[0]&&p<=D.gap[1];}
 function keyMap(rows,mx,posf){var o={};(rows||[]).forEach(function(r){if(r.matrix!==mx)return;o[r.time+"|"+posf(r)]=r;});return o;}
 function covSet(mx){var s={};(D.covered||[]).forEach(function(r){if(r.matrix===mx)s[r.time+"|"+r.pos]=1;});return s;}
 function esc(s){return String(s).replace(/</g,"&lt;");}
 function h(t,a,inner){var s="<"+t; for(var k in a)s+=" "+k+"=\\\""+a[k]+"\\\""; return s+">"+(inner||"")+"</"+t.split(" ")[0]+">";}

 function renderGrid(mx,host){
   var bm=keyMap(D.base,mx,function(r){return r.pos;}), cov=covSet(mx), out="";
   var head="<tr><th class=rl></th>"; for(var p=1;p<=N;p++)head+="<th class=ph>"+p+"</th>"; head+="<th class=rl></th></tr>";
   var seq="<tr><th class=rl style=\\"color:var(--mut)\\">5\\u2032\\u21923\\u2032</th>";
   for(var p=1;p<=N;p++)seq+="<td class=\\"sc"+(inGap(p)?" g":"")+"\\">"+esc(D.seq[p-1])+"</td>";
   seq+="<th class=rl style=\\"color:var(--mut);font-size:9.5px\\">parent</th></tr>";
   var body="";
   D.times.forEach(function(tm){
     var row="<tr><th class=rl>"+esc(tm)+"</th>";
     for(var p=1;p<=N;p++){
       var cell=bm[tm+"|"+p];
       if(cell){var cl=CLS[cell.cls],col=cvar(cl[0]),op=alpha(cell.pct,7);
         row+="<td class=c data-tip=\\""+tipB(cell,p,tm,mx)+"\\"><div class=bx style=\\"background:"+rgba(col,op)+"\\">"+cl[1]+"</div>"+(cell.ms2?"<span class=d></span>":"")+"</td>";
       } else if(cov[tm+"|"+p]){row+="<td class=\\"c cov\\" data-tip=\\""+tipCov(p,tm,mx)+"\\"><div class=bx>\\u00b7</div></td>";}
       else row+="<td class=c><div class=bx></div></td>";
     }
     var pr=parentPct(mx,tm);
     row+="<th class=rl>"+(pr==null?"":"<span class=trk><span style=\\"width:"+Math.max(0,Math.min(100,pr)).toFixed(0)+"%\\"></span></span><span class=p>"+pr.toFixed(1)+"%</span>")+"</th></tr>";
     body+=row;
   });
   host.innerHTML="<div class=scrollx><table>"+head+seq+body+"</table></div>";
   host.querySelectorAll("[data-tip]").forEach(function(td){bind(td,td.getAttribute("data-tip"));td.removeAttribute("data-tip");});
 }
 function renderCleave(mx,host){
   var cm=keyMap(D.cleave,mx,function(r){return r.link;}), any=(D.cleave||[]).some(function(r){return r.matrix===mx;});
   if(!any){host.innerHTML="<div class=empty>No backbone-cleavage products above threshold in this matrix.</div>";return;}
   var head="<tr><th class=rl></th>"; for(var l=1;l<N;l++)head+="<th class=ph>"+l+"</th>"; head+="</tr>";
   var body="";
   D.times.forEach(function(tm){var row="<tr><th class=rl>"+esc(tm)+"</th>";
     for(var l=1;l<N;l++){var c=cm[tm+"|"+l];
       if(c){var cl=CLS[c.cls],op=alpha(c.pct,6);
         row+="<td class=c data-tip=\\""+tipC(c,l,tm,mx)+"\\"><div class=bx style=\\"background:"+rgba(cvar(cl[0]),op)+"\\">"+cl[1]+"</div>"+(c.ms2?"<span class=d></span>":"")+(c.cyc?"<span class=cp>&gt;p</span>":"")+"</td>";
       } else row+="<td class=c><div class=bx></div></td>";
     } row+="</tr>"; body+=row;});
   host.innerHTML="<div class=scrollx><table>"+head+body+"</table></div>";
   host.querySelectorAll("[data-tip]").forEach(function(td){bind(td,td.getAttribute("data-tip"));td.removeAttribute("data-tip");});
 }
 function renderConj(mx,host){
   var rows=(D.conj||[]).filter(function(r){return r.matrix===mx;});
   if(!rows.length){host.innerHTML="<div class=empty>No conjugate catabolites (parent is not GalNAc-conjugated, or none detected).</div>";return;}
   var labels=[]; rows.forEach(function(r){if(labels.indexOf(r.label)<0)labels.push(r.label);});
   var m={}; rows.forEach(function(r){m[r.time+"|"+r.label]=r;});
   var head="<tr><th class=rl></th>"; labels.forEach(function(l){head+="<th class=ph style=\\"width:92px;font-size:9.5px\\">"+esc(l.replace(/^.*?(-\\d GalNAc|scaffold|de-?conjug.*|linker.*)$/i,"$1"))+"</th>";}); head+="</tr>";
   var body="";
   D.times.forEach(function(tm){var row="<tr><th class=rl>"+esc(tm)+"</th>";
     labels.forEach(function(l){var c=m[tm+"|"+l];
       if(c){var op=alpha(c.pct,100);
         row+="<td class=c style=width:92px data-tip=\\""+tipGn(c,tm,mx)+"\\"><div class=bx style=\\"background:"+rgba(cvar("--gn"),op)+";inset:2px 8px\\">"+c.pct.toFixed(0)+"%</div></td>";
       } else row+="<td class=c style=width:92px><div class=bx></div></td>";
     }); row+="</tr>"; body+=row;});
   host.innerHTML="<div class=scrollx><table>"+head+body+"</table></div>";
   host.querySelectorAll("[data-tip]").forEach(function(td){bind(td,td.getAttribute("data-tip"));td.removeAttribute("data-tip");});
 }
 function parentPct(mx,tm){var r=(D.parent||[]).find(function(x){return x.matrix===mx&&x.time===tm;});return r?r.pct:null;}
 function tipB(c,p,tm,mx){var cl=CLS[c.cls];return "<b>"+cl[2]+"</b><br><span class=k>residue</span> "+esc(D.seq[p-1])+p+(inGap(p)?" (DNA gap)":" (wing)")+"<hr><span class=k>species</span> "+esc(c.species)+"<br><span class=k>abundance</span> "+c.pct.toFixed(2)+"% of parent<br>"+(isFinite(c.ppm)?"<span class=k>mass error</span> "+(+c.ppm).toFixed(1)+" ppm<br>":"")+"<span class=k>MS2</span> "+(c.ms2?"confirmed":"MS1 only")+"<br><span class=k>"+esc(tm)+" \\u00b7 "+esc(mx)+"</span>";}
 function tipC(c,l,tm,mx){var cl=CLS[c.cls];return "<b>"+cl[2]+"</b> <span class=k>@ linkage "+l+"</span><hr><span class=k>product</span> "+esc(c.species)+(c.cyc?" (2\\u2032,3\\u2032-cyclic phosphate)":"")+"<br><span class=k>abundance</span> "+c.pct.toFixed(2)+"% of parent<br><span class=k>MS2</span> "+(c.ms2?"confirmed":"MS1 only")+"<br><span class=k>"+esc(tm)+" \\u00b7 "+esc(mx)+"</span>";}
 function tipGn(c,tm,mx){return "<b>"+esc(c.label)+"</b><hr><span class=k>abundance</span> "+c.pct.toFixed(1)+"% of parent<br><span class=k>MS2</span> "+(c.ms2?"confirmed":"MS1 only")+"<br><span class=k>"+esc(tm)+" \\u00b7 "+esc(mx)+"</span>";}
 function tipCov(p,tm,mx){return "<b>pos "+p+" \\u00b7 "+esc(D.seq[p-1])+"</b> <span class=k>("+(inGap(p)?"DNA gap":"wing")+")</span><hr><span class=k>covered, no modification</span><br>"+esc(tm)+" \\u00b7 "+esc(mx);}

 function renderAll(mx){
   renderGrid(mx,root.querySelector("#g_"+root.id));
   renderCleave(mx,root.querySelector("#c_"+root.id));
   renderConj(mx,root.querySelector("#j_"+root.id));
 }
 var tabs=root.querySelectorAll(".tab");
 tabs.forEach(function(b){b.addEventListener("click",function(){tabs.forEach(function(x){x.setAttribute("aria-selected",x===b);});renderAll(b.dataset.m);});});
 renderAll(D.matrices[0]);
 if(window.matchMedia){var mq=matchMedia("(prefers-color-scheme:dark)"); if(mq.addEventListener)mq.addEventListener("change",function(){var sel=root.querySelector(".tab[aria-selected=true]");renderAll(sel?sel.dataset.m:D.matrices[0]);});}
})();
</script>'

#' Render the degradation position map as a self-contained HTML widget
#'
#' Produces an HTML fragment (scoped under `.dmap`) that draws the model
#' from [degradation_position_model()] as three aligned panels (base
#' modifications, backbone cleavage, GalNAc conjugate catabolism) with a
#' matrix selector, a parent-remaining track, hover tooltips, and a legend.
#' Self-contained: no external scripts or styles, safe to embed in the Shiny
#' app (via `shiny::HTML`) or the study report.
#'
#' @param model A [degradation_position_model()] result.
#' @param id A unique element id (so several widgets can coexist).
#' @return A single HTML string. An empty-state message when the model has
#'   no placeable cells.
#' @seealso [degradation_position_model()]
#' @export
degradation_map_html <- function(model, id = "dmap1") {
  mats <- model$matrices
  has_cells <- (nrow(model$base) + nrow(model$cleave) + nrow(model$conj)) > 0
  if (length(mats) == 0 || !has_cells) {
    return(paste0('<div class="dmap"><div class="empty" style="padding:18px">',
                  'No positioned degradation products yet. Run batch processing (or import a ',
                  'data matrix) with sample metadata so metabolites can be placed by position, ',
                  'timepoint, and matrix.</div></div>'))
  }
  payload <- list(
    oligo = model$oligo, seq = as.character(model$seq),
    gap = if (is.null(model$gap)) NULL else as.integer(model$gap),
    times = as.character(model$times), matrices = as.character(mats),
    base = model$base, cleave = model$cleave, conj = model$conj,
    covered = model$covered, parent = model$parent)
  json <- as.character(jsonlite::toJSON(payload, auto_unbox = TRUE,
                                        dataframe = "rows", na = "null", null = "null"))
  tabs <- paste0(vapply(seq_along(mats), function(i) sprintf(
    '<button class="tab" role="tab" data-m="%s" aria-selected="%s">%s</button>',
    mats[i], tolower(i == 1), mats[i]), character(1)), collapse = "")
  legend <- paste0(
    '<div class="panel"><h3>Legend</h3><div class="legend">',
    '<div class="lg"><span class="t">Base modification</span>',
    '<span class="lr"><span class="sw" style="background:var(--ox)">ox</span> Oxidation</span>',
    '<span class="lr"><span class="sw" style="background:var(--aB)">aB</span> Depurination / abasic</span>',
    '<span class="lr"><span class="sw" style="background:var(--dA)">dA</span> Deamination</span></div>',
    '<div class="lg"><span class="t">Backbone cleavage</span>',
    '<span class="lr"><span class="sw" style="background:var(--e3)">3&prime;</span> 3&prime;-exonuclease</span>',
    '<span class="lr"><span class="sw" style="background:var(--e5)">5&prime;</span> 5&prime;-exonuclease</span>',
    '<span class="lr"><span class="sw" style="background:var(--en)">en</span> Endonuclease</span></div>',
    '<div class="lg"><span class="t">Conjugate / coverage</span>',
    '<span class="lr"><span class="sw" style="background:var(--gn)">Gn</span> GalNAc loss</span>',
    '<span class="lr"><span class="sw cov">&nbsp;</span> Covered, no modification</span></div>',
    '<div class="lg"><span class="t">Confidence</span>',
    '<span class="lr"><span class="mkdot"></span>&nbsp; MS2-confirmed</span>',
    '<span class="lr" style="color:var(--mut)">Opacity &prop; % of parent &middot; &gt;p = cyclic phosphate</span></div>',
    '</div></div>')
  body <- paste0(
    '<div class="dmap" id="', id, '">',
    '<div class="tabs" role="tablist">', tabs, '</div>',
    '<div class="panel"><h3>In-place base modifications &middot; residue &times; timepoint</h3><div id="g_', id, '"></div></div>',
    '<div class="panel"><h3>Backbone cleavage &middot; linkage &times; timepoint</h3><div id="c_', id, '"></div></div>',
    '<div class="panel"><h3>GalNAc conjugate catabolism &middot; state &times; timepoint</h3><div id="j_', id, '"></div></div>',
    legend,
    '<div class="tip"></div></div>')
  js <- sub("__DATA__", json, .DM_WIDGET_JS, fixed = TRUE)
  js <- sub("__ID__", id, js, fixed = TRUE)
  paste0(.DM_WIDGET_CSS, body, js)
}
