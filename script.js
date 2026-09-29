// ---- Panels Toggle ----

const panels = document.querySelectorAll(".panel")
const navLinks = document.querySelectorAll("#nav a")

function showPanel(id) {
  panels.forEach((p) => p.classList.toggle("active", p.id === id))
  navLinks.forEach((a) =>
    a.classList.toggle("active", a.getAttribute("href") === `#${id}`),
  )
}

navLinks.forEach((link) => {
  link.addEventListener("click", (e) => {
    const href = link.getAttribute("href")
    if (href.startsWith("#") && href.length > 1) {
      e.preventDefault()
      const id = href.slice(1)
      showPanel(id)
      history.pushState(null, "", href)
    }
  })
})

function initFromHash() {
  const id = location.hash ? location.hash.slice(1) : "home"
  showPanel(document.getElementById(id) ? id : "home")
}

window.addEventListener("popstate", initFromHash)
initFromHash()

//----- contact copy email------

async function copyEmail() {
  const email = "asandatope@gmail.com"
  const message = document.querySelector("#copy-message")

  try {
    await navigator.clipboard.writeText(email)
    message.textContent = "Email copied!"

    setTimeout(() => {
      message.textContent = ""
    }, 3000) // disappears after 3 seconds
  } catch (error) {
    console.log("ERR", error)
    message.textContent = "Couldn't copy email address"

    setTimeout(() => {
      message.textContent = ""
    }, 3000)
  }
}

const emailButton = document.querySelector("#email-button")
emailButton.addEventListener("click", copyEmail)

// ---- Footer year ----
document.getElementById("year").textContent = new Date().getFullYear()

// ---- Load and render projects from projects.json ----
async function loadProjects() {
  const tbody = document.querySelector("#project-list tbody")
  if (!tbody) return

  try {
    const res = await fetch("projects.json")
    if (!res.ok) throw new Error(`Failed to load projects.json (${res.status})`)
    const projects = await res.json()

    if (!projects.length) {
      tbody.innerHTML = `<tr><td colspan="6" class="text-sm text-ink/50 dark:text-panelText/50">No projects yet — check back soon.</td></tr>`
      return
    }

    projects.sort((a, b) => (b.year || 0) - (a.year || 0))
    tbody.innerHTML = projects.map(renderProjectRow).join("")
  } catch (err) {
    console.error(err)
    tbody.innerHTML = `<tr><td colspan="6" class="text-sm text-red-600 dark:text-red-400">Couldn't load projects right now.</td></tr>`
  }
}

function renderProjectRow(p) {
  const linkHtml = p.link
    ? `<a href="${escapeHtml(p.link)}" class="btn rounded-none" data-variant="outline" data-size="sm" target="_blank" rel="noopener">
                   <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-arrow-up-right preview-icon"><path d="M7 7h10v10"/><path d="M7 17 17 7"/></svg></a>`
    : `<span class="text-xs  text-ink/40 dark:text-panelText/40">—</span>`

  const techHtml = p.tech?.length
    ? `<div class="text-xs text-ink/50 dark:text-panelText/50 mt-1">${escapeHtml(p.tech.join(", "))}</div>`
    : ""

  return `
    <tr>
      <td class="font-medium text-wrap">${escapeHtml(p.title)}</td>
      <td class="text-wrap">
        ${escapeHtml(p.description || "")}
        ${techHtml}
      </td>
      <td class="truncate">${escapeHtml(p.employer || "")}</td>
      <td>${escapeHtml(p.status || "")}</td>
      <td>${escapeHtml(String(p.year || ""))}</td>
      <td class="text-end">${linkHtml}</td>
    </tr>
  `
}

function escapeHtml(str) {
  const div = document.createElement("div")
  div.textContent = str
  return div.innerHTML
}

loadProjects()
