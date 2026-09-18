(() => {
  "use strict";


  /* ================================================================
     STATE
     ================================================================ */

  let activeNode = null;

  let currentTree = null;
  let currentPanel = null;
  let currentTitle = null;
  let currentBody = null;
  let currentLink = null;

  let birdTimer = null;

  /*
   * Keeps track of the nodes that have already received listeners.
   * This prevents duplicate listeners if Quartz reinitializes the page.
   */
  const boundNodes = new WeakSet();


  /* ================================================================
     PANEL DATA
     ================================================================ */

  const panelData = {

    personal: {
      title: "Personal Information",
      body:
        "The personal and academic foundation of the researcher: " +
        "contact information, education, academic records, and " +
        "professional presence.",
      label: "Open personal information →",
      url: ""
    },

    CV: {
      title: "CV",
      label: "Curriculum Vitae",
      body: 
"This latin word translates to 'course of life'. Also, a fun fact is that a letter from Leonardo Da Vinci to the Duke of Milan in 1498 is considered to be one of the earliest records resembling a modern CV (but it read more like a cover letter).",  
      url: ""
    },

    grades: {
      title: "Academic Records",
      body:
        "Academic transcripts, grades, qualifications, and " +
        "other educational records.",
      label: "Open academic records →",
      url: ""
    },

    motivation: {
      title: "Intellectual Motivation",
      body:
        "The questions and observations that motivate the research " +
        "programme and connect otherwise different biological systems.",
      label: "Research philosophy →",
      url: ""
    },

    "motivation-detail": {
      title: "Research Motivation",
      body:
        "Understanding how physical context influences biological " +
        "function across scales, from molecular interactions to " +
        "cells, tissues, and disease progression.",
      label: "Read motivation →",
      url: ""
    },

    principles: {
      title: "Guiding Principles",
      body:
        "Different biological systems often solve different problems " +
        "using similar physical and mechanistic principles.",
      label: "Read guiding principles →",
      url: ""
    },



    mechanobiology: {
      title: "Mechanobiology",
      body:
        "An introduction to how mechanical forces, stiffness, " +
        "geometry, and physical microenvironments influence " +
        "biological signalling, organisation, behaviour, and fate.",
      label: "Open mechanobiology section →",
      url: "https://asharrais90.github.io/Research-Projects-Ashar-Rais-/"
    },

    project1: {
      title: "Project 1 — Thymic Selection",
      body:
        "Explore how mechanical and biochemical properties of " +
        "the thymic environment may influence thymocyte selection " +
        "and T-cell development.",
      label: "Open project →",
      url: "https://asharrais90.github.io/Research-Projects-Ashar-Rais-/thymic-selection-mechanobiology/"
    },

    project2: {
      title: "Project 2 — Pre-metastatic Niche",
      body:
        "Investigate how physical and microenvironmental changes " +
        "can contribute to the establishment of a pre-metastatic niche.",
      label: "Open project →",
      url: ""
    },

    project3: {
      title: "Project 3 — Computational Biology",
      body:
        "Mathematical modelling, simulation, quantitative analysis, " +
        "and reproducible computational workflows for studying " +
        "biological systems.",
      label: "Open project →",
      url: ""
    },

  };

const ROOT_CONFIGS = [
  {
    startX: 540,
    endX: 320,
    depth: 155,
    thickness: 36,
    curve: -30,
    seed: 41,
    label: "Guiding Principles",
    panel: "principles"
  },
  {
    startX: 560,
    endX: 600,
    depth: 180,
    thickness: 43,
    curve: 5,
    seed: 47,
    label: "CV",
    panel: "CV"
  },
  {
    startX: 585,
    endX: 850,
    depth: 145,
    thickness: 35,
    curve: 30,
    seed: 62,
    label: "Intellectual Motivation & Background",
    panel: "personal"
  }
];


  /* ================================================================
     PROCEDURAL FOLIAGE
     ================================================================ */

  function generateFoliage(canopy) {

    if (!canopy) {
      return;
    }

    /*
     * Clear any foliage generated during a previous initialization.
     */

    canopy.innerHTML = "";


    /*
     * Each object controls one foliage mass.
     *
     * x, y       = position
     * rx, ry     = overall size
     * roughness  = local irregularity
     * seed       = deterministic shape variation
     */

    const clusters = [

      {x: 190, y: 230, rx: 185, ry: 135, roughness: 0.45, seed: 11},

      {x: 350, y: 175, rx: 200, ry: 165, roughness: 0.49, seed: 27},

      {x: 550, y: 155, rx: 210, ry: 160, roughness: 0.48, seed: 49},

      {x: 750, y: 178, rx: 210, ry: 165, roughness: 0.50, seed: 63},

      {x: 915, y: 235, rx: 185, ry: 165, roughness: 0.47, seed: 81},

      {x: 300, y: 290, rx: 230, ry: 105, roughness: 0.50, seed: 101},

      {x: 530, y: 295, rx: 250, ry: 115, roughness: 0.47, seed: 121},

      {x: 770, y: 295, rx: 230, ry: 105, roughness: 0.49, seed: 147}

    ];


    // 1. Create ONE group to hold ALL the tufts
    const foliageGroup = document.createElementNS("http://www.w3.org/2000/svg",
      "g");
    foliageGroup.setAttribute("id", "canopyGroup");

    // 2. Loop through the clusters (EXACTLY as before)
    clusters.forEach((config) => {

      const path = document.createElementNS("http://www.w3.org/2000/svg", "path");

      path.setAttribute("class", "foliage-cluster"
      );

      path.setAttribute(
        "d",
        createOrganicFoliagePath(config)
      );

      // 3. Append the path to the GROUP, not directly to the canopy
      foliageGroup.appendChild(path);

    });

    // 4. Append the complete group to the canopy
      canopy.appendChild(foliageGroup);
  }


  /* ================================================================
     ORGANIC FOLIAGE PATH
     ================================================================ */

  function createOrganicFoliagePath(config) {

    const {x, y, rx, ry, roughness, seed
    } = config;


    const points = [];

    /*
     * More points = finer canopy boundary.
     */

    const pointCount = 80;


    for (let i = 0; i < pointCount; i++) {

      const angle =
        (Math.PI * 2 * i) /
        pointCount;


      const noise =
        deterministicNoise(i, seed);


      /*
       * Large-scale shape.
       */

      const largeShape =
        0.10 * Math.sin(angle * 2 + seed * 0.11);


      /*
       * Medium-scale foliage lobes.
       */

      const mediumLobes =
        0.095 * Math.sin(angle * 7 + seed * 0.31);


      /*
       * Small-scale edge irregularity.
       *
       * roughness remains available per cluster.
       */

      const smallTexture = (noise - 0.55) * roughness * 0.45;


      const variation = 1 + largeShape + mediumLobes + smallTexture;


      const px = x + Math.cos(angle) * rx * variation;


      const py = y + Math.sin(angle) * ry * variation;


      points.push({x: px, y: py});

    }


    /*
     * Build a smooth closed path through the irregular points.
     */

    let path =
      `M ${points[0].x} ${points[0].y}`;


    for (let i = 0; i < pointCount; i++) {

      const current =
        points[i];

      const next =
        points[(i + 1) % pointCount];


      const midpointX =
        (current.x + next.x) / 2;

      const midpointY =
        (current.y + next.y) / 2;


      path +=
        ` Q ${current.x} ${current.y}` +
        ` ${midpointX} ${midpointY}`;

    }


    path += " Z";

    return path;
  }

  /* ================================================================
     ORGANIC ROOT PATH & GROUND
     ================================================================ */


function generateRoots(rootGroup) {
  if (!rootGroup) return;

  rootGroup.innerHTML = "";

  ROOT_CONFIGS.forEach((config) => {
    const path = createOrganicRootPath(config);
    path.classList.add("root-shape");
    rootGroup.appendChild(path);
  });
}


function generateRootLabels(labelGroup) {
  if (!labelGroup) return;

  labelGroup.innerHTML = "";

  ROOT_CONFIGS.forEach((config) => {
    const tipX = config.endX;
    const tipY = 790 + config.depth + 28;   // 28px below the tip

    const group = document.createElementNS(
      "http://www.w3.org/2000/svg",
      "g"
    );
    group.classList.add("root-label");
    group.setAttribute("tabindex", "0");
    group.setAttribute("role", "button");
    group.dataset.panel = config.panel;

    const text = document.createElementNS(
      "http://www.w3.org/2000/svg",
      "text"
    );
    text.classList.add("root-label-text");
    text.setAttribute("x", tipX);
    text.setAttribute("y", tipY);
    text.setAttribute("text-anchor", "middle");
    text.textContent = config.label;

    group.appendChild(text);
    labelGroup.appendChild(group);
  });
}


function createOrganicRootPath(config) {
  const {
    startX,
    endX,
    depth,
    thickness,
    curve,
    seed
  } = config;

  const startY = 790;
  const dx = endX - startX;

  const control1X = startX + dx * 0.18;
  const control1Y = startY + depth * 0.18;

  const control2X = startX + dx * 0.55;
  const control2Y = startY + depth * 0.72 + curve;

  const endY = startY + depth;

  const angle = Math.atan2(endY - startY, endX - startX);

  const normalX = Math.cos(angle + Math.PI / 2);
  const normalY = Math.sin(angle + Math.PI / 2);

  const baseOffsetX = normalX * thickness * 0.5;
  const baseOffsetY = normalY * thickness * 0.5;

  const tipWidth = thickness * 0.08;

  const tipOffsetX = normalX * tipWidth;
  const tipOffsetY = normalY * tipWidth;

  const bend = deterministicNoise(seed, 0) * 18 - 9;

  const d = `
    M ${startX - baseOffsetX}
      ${startY - baseOffsetY}

    C ${control1X - baseOffsetX}
      ${control1Y - baseOffsetY},

      ${control2X + bend - baseOffsetX * 0.5}
      ${control2Y - baseOffsetY * 0.5},

      ${endX - tipOffsetX}
      ${endY - tipOffsetY}

    C ${endX + tipOffsetX}
      ${endY + tipOffsetY},

      ${control2X + bend + baseOffsetX * 0.5}
      ${control2Y + baseOffsetY * 0.5},

      ${control1X + baseOffsetX}
      ${control1Y + baseOffsetY}

    C ${startX + baseOffsetX}
      ${startY + baseOffsetY},

      ${startX}
      ${startY + thickness * 0.45},

      ${startX - baseOffsetX}
      ${startY - baseOffsetY}

    Z
  `;

  const path = document.createElementNS(
    "http://www.w3.org/2000/svg",
    "path"
  );

  path.setAttribute("d", d);

  return path;
}


function generateGroundTexture(groundGroup) {
  if (!groundGroup) return;

  // Remove any previously generated texture
  groundGroup.querySelectorAll(".ground-tuft").forEach((n) => n.remove());

  const tufts = [];

  // Spread tufts along the visible ground line (x = 160 to 1040)
  for (let i = 0; i < 42; i++) {
    const seed = i * 13 + 7;
    const x = 180 + deterministicNoise(i, 91) * 860;
    const y = 800 + deterministicNoise(i, 137) * 6 - 2;
    const height = 4 + deterministicNoise(i, 173) * 10;
    const width = 1.5 + deterministicNoise(i, 211) * 1.5;
    const lean = (deterministicNoise(i, 251) - 0.5) * 6;

    tufts.push({ x, y, height, width, lean, seed });
  }

  tufts.forEach((t) => {
    const path = document.createElementNS(
      "http://www.w3.org/2000/svg",
      "path"
    );

    const { x, y, height, width, lean } = t;

    path.setAttribute(
      "d",
      `
        M ${x - width} ${y}
        Q ${x + lean * 0.3} ${y - height * 0.6},
          ${x + lean} ${y - height}
        Q ${x + width * 0.4} ${y - height * 0.4},
          ${x + width} ${y}
        Z
      `
    );

    path.classList.add("ground-tuft");
    groundGroup.appendChild(path);
  });
}
 

 /* ================================================================
     COLOR TRANSITION ANIMATION
     ================================================================ */

let colorCycleFrame = null;

function startTreeGradientAnimation() {
  if (colorCycleFrame !== null) return;
  if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) return;

  let startTime = null;

  function cycleColors(timestamp) {
    if (!startTime) startTime = timestamp;
    
    // Gentle breathing speed: Full cycle (brown → green → brown) takes ~2 minutes.
    // Change the 0.015 to adjust speed. Larger = faster.
    const elapsed = (timestamp - startTime) / 1000;
    const progress = (Math.sin(elapsed * 0.255) + 1) / 2; // Progress 0 to 1 and back

    // ----- FORESTY PALETTE -----
    // Hue: 35 (Earthy Brown) to 115 (Deep Forest Green)
    const hue = 35 + (progress * 80);
    
    // Saturation: 40% (muted, academic) to 65% (richer, alive)
    const saturation = 40 + (progress * 15);
    
    // Lightness: 30% (dark, grounding) to 45% (brighter, growing)
    const lightness = 30 + (progress * 15);

    // Create the final color string
    const currentColor = `hsl(${hue}, ${saturation}%, ${lightness}%)`;

    const trunkGroup = document.getElementById("trunkGroup");
    if (trunkGroup) {
      trunkGroup.setAttribute('fill', currentColor);
    }

    const canopyGroup = document.getElementById("canopyGroup");
    if (canopyGroup) {
      canopyGroup.setAttribute('fill', currentColor);
    }

    const groundGroup = document.getElementById("groundGroup");
    if (groundGroup) {
       groundGroup.setAttribute("fill", currentColor);
    }

    // ----- 3. (OPTIONAL) If you have a separate gradient for highlights, keep it static -----
    // You can delete the old `treeColourGradient` stops code entirely now.

    colorCycleFrame = requestAnimationFrame(cycleColors);
  }

  colorCycleFrame = requestAnimationFrame(cycleColors);
}


 /* ================================================================
     DETERMINISTIC NOISE
     ================================================================ */

  function deterministicNoise(index, seed) {

    const value =
      Math.sin(index * 12.9898 + seed * 78.233) * 43758.5453;


    return (value - Math.floor(value));

  }


  /* ================================================================
     OPEN INFORMATION PANEL
     ================================================================ */

  function openPanel(panelId, node) {

    const data =
      panelData[panelId];


    if (
      !data ||
      !currentPanel ||
      !currentTitle ||
      !currentBody ||
      !currentLink
    ) {
      return;
    }


    currentTitle.textContent =
      data.title;

    currentBody.textContent =
      data.body;


    if (data.url) {

      currentLink.hidden = false;

      currentLink.href =
        data.url;

      currentLink.textContent =
        data.label;


      if (
        /^https?:\/\//i.test(
          data.url
        )
      ) {

        currentLink.target =
          "_blank";

        currentLink.rel =
          "noopener noreferrer";

      } else {

        currentLink.removeAttribute(
          "target"
        );

        currentLink.removeAttribute(
          "rel"
        );

      }

    } else {

      currentLink.hidden = true;

      currentLink.removeAttribute(
        "href"
      );

      currentLink.textContent =
        "";

    }


    currentPanel.classList.add(
      "is-open"
    );

    currentPanel.setAttribute(
      "aria-hidden",
      "false"
    );


    if (
      activeNode &&
      activeNode !== node
    ) {

      activeNode.setAttribute(
        "aria-expanded",
        "false"
      );

    }


    activeNode = node;

    node.setAttribute(
      "aria-expanded",
      "true"
    );

  }


  /* ================================================================
     CLOSE INFORMATION PANEL
     ================================================================ */

  function closePanel() {

    if (!currentPanel) {
      return;
    }


    currentPanel.classList.remove(
      "is-open"
    );

    currentPanel.setAttribute(
      "aria-hidden",
      "true"
    );


    if (activeNode) {

      activeNode.setAttribute(
        "aria-expanded",
        "false"
      );

      /*
       * Do not force focus back if the element was removed by
       * Quartz navigation.
       */

      if (
        document.contains(activeNode)
      ) {

        activeNode.focus({
          preventScroll: true
        });

      }

      activeNode = null;

    }

  }


  /* ================================================================
     INITIALIZE / RE-INITIALIZE TREE
     ================================================================ */

  function initializeResearchTree() {

    const tree =
      document.getElementById(
        "research-tree"
      );

    const panel =
      document.getElementById(
        "info-panel"
      );

    const title =
      document.getElementById(
        "panel-title"
      );

    const body =
      document.getElementById(
        "panel-body"
      );

    const link =
      document.getElementById(
        "panel-link"
      );

    const closeButton =
      document.querySelector(
        ".panel-close"
      );

    const canopy =
      document.getElementById(
        "canopy"
      );

    const roots =
      document.getElementById(
        "roots"
      );
    
    const groundGroup = document.getElementById("groundGroup");
    
    const rootLabels = tree.querySelector("#root-labels");

    const ground = document.getElementById("ground");
    
    if (canopy) {
      generateFoliage(canopy);
    }

    if (roots) {
      generateRoots(roots);
      if (rootLabels) generateRootLabels(rootLabels);
    }

    if (ground) {
      generateGroundTexture(ground);
    }
    

    /*
     * This simply means that the current Quartz page is not
     * the Central Hub. That is normal.
     */

    if (
      !tree ||
      !panel ||
      !title ||
      !body ||
      !link ||
      !closeButton ||
      !canopy ||
      !roots 
    ) {

      currentTree = null;
      currentPanel = null;
      currentTitle = null;
      currentBody = null;
      currentLink = null;
      activeNode = null;

      return;
    }


    currentTree =
      tree;

    currentPanel =
      panel;

    currentTitle =
      title;

    currentBody =
      body;

    currentLink =
      link;


    /*
     * Generate the foliage.
     */

    generateFoliage(
      canopy
    );


    /*
     * Reset the current interaction state.
     */

    activeNode = null;


    /*
     * Bind all interactive tree nodes.
     *
     * WeakSet prevents duplicate listeners on the same DOM node.
     */

    tree
      .querySelectorAll("[data-panel], [data-url]")
      .forEach((node) => {

        if (
          boundNodes.has(node)
        ) {
          return;
        }


        boundNodes.add(node);


        node.addEventListener(
          "click",
          (event) => {

            event.stopPropagation();

              if (node.dataset.url) {
    window.open(
      node.dataset.url,
      "_blank",
      "noopener,noreferrer"
    );
    return;
  }  

            const panelId =
              node.dataset.panel;


            if (
              activeNode === node
            ) {

              closePanel();

              return;
            }


            openPanel(
              panelId,
              node
            );

          }
        );


        node.addEventListener(
          "keydown",
          (event) => {

            if (
              event.key === "Enter" ||
              event.key === " "
            ) {

              event.preventDefault();

              node.click();

            }
      });

if (node.dataset.tooltip) {
  const tooltip = document.getElementById("canopy-tooltip");
  const stage = document.getElementById("research-tree");

  const showTooltip = () => {
    if (!tooltip || !stage) return;

    // Position the tooltip under the node
    const stageRect = stage.getBoundingClientRect();
    const nodeRect = node.getBoundingClientRect();

    tooltip.style.left =
      `${nodeRect.left - stageRect.left + nodeRect.width / 2}px`;
    tooltip.style.top =
      `${nodeRect.bottom - stageRect.top + 10}px`;

    tooltip.classList.add("is-visible");
    tooltip.setAttribute("aria-hidden", "false");
  };

  const hideTooltip = () => {
    if (!tooltip) return;
    tooltip.classList.remove("is-visible");
    tooltip.setAttribute("aria-hidden", "true");
  };

  node.addEventListener("mouseenter", showTooltip);
  node.addEventListener("mouseleave", hideTooltip);
  node.addEventListener("focus", showTooltip);
  node.addEventListener("blur", hideTooltip);
}
});

    /*
     * Bind this page's close button.
     */

    if (
      !closeButton.dataset.bound
    ) {

      closeButton.dataset.bound =
        "true";

      closeButton.addEventListener(
        "click",
        (event) => {

          event.stopPropagation();

          closePanel();

        }
      );

    }

  }


  /* ================================================================
     BIRDS
     ================================================================ */

  function createBird() {

    if (!currentTree) {
      return;
    }


    const bird =
      document.createElement(
        "div"
      );


    bird.className =
      "bird bird-extra";

    bird.setAttribute(
      "aria-hidden",
      "true"
    );


    bird.innerHTML = `
      <svg viewBox="0 0 120 60">
        <path
          d="
            M8 34
            Q30 8 58 28
            Q86 8 112 34
            Q84 20 58 42
            Q31 20 8 34Z
          "
        />
      </svg>
    `;


    bird.style.top =
      `${8 + Math.random() * 25}%`;


    bird.style.width =
      `${18 + Math.random() * 18}px`;


    bird.style.height =
      `${10 + Math.random() * 10}px`;


    bird.style.opacity =
      `${0.12 + Math.random() * 0.22}`;


    bird.style.animationDuration =
      `${24 + Math.random() * 18}s`;


    currentTree.appendChild(
      bird
    );


    window.setTimeout(
      () => {

        if (
          bird.isConnected
        ) {

          bird.remove();

        }

      },
      48000
    );

  }


  function startBirdTimer() {

    /*
     * Only create ONE interval.
     */

    if (birdTimer !== null) {
      return;
    }


    birdTimer =
      window.setInterval(
        () => {

          /*
           * Don't generate birds if we're no longer on the
           * Central Hub.
           */

          if (
            !currentTree
          ) {
            return;
          }


          if (
            window.matchMedia(
              "(prefers-reduced-motion: reduce)"
            ).matches
          ) {

            return;

          }


          createBird();

        },
        19000
      );

  }


  /* ================================================================
     GLOBAL KEYBOARD HANDLER
     ================================================================ */

  document.addEventListener(
    "keydown",
    (event) => {

      if (
        event.key === "Escape" &&
        currentPanel &&
        currentPanel.classList.contains(
          "is-open"
        )
      ) {

        closePanel();

      }

    }
  );


  /* ================================================================
     GLOBAL OUTSIDE-CLICK HANDLER
     ================================================================ */

  document.addEventListener(
    "click",
    (event) => {

      if (
        !currentPanel ||
        !currentPanel.classList.contains(
          "is-open"
        )
      ) {

        return;

      }


      if (
        !currentPanel.contains(
          event.target
        ) &&
        currentTree &&
        !currentTree.contains(
          event.target
        )
      ) {

        closePanel();

      }

    }
  );


  /* ================================================================
     INITIAL PAGE LOAD
     ================================================================ */

  if (
    document.readyState ===
    "loading"
  ) {

    document.addEventListener(
      "DOMContentLoaded",
      () => {
        initializeResearchTree();
        startBirdTimer();
        startTreeGradientAnimation();
      }, { once: true });

  } else {
    initializeResearchTree();
    startBirdTimer();
    startTreeGradientAnimation();
  }


  /* ================================================================
     QUARTZ SPA NAVIGATION
     ================================================================ */

  document.addEventListener(
    "nav",
    () => {

      initializeResearchTree();
      startBirdTimer();
      startTreeGradientAnimation();
    }
  );

})();